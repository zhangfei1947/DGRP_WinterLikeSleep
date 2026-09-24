#!/usr/bin/env Rscript
# Rscript scripts/extract_profile_phenotypes.R [data_root] [output_root] [batch_regex]
suppressPackageStartupMessages({library(data.table);library(damr);library(behavr)})
source("scripts/profile_phenotype_helpers.R")
setDTthreads(2)
args <- commandArgs(TRUE)
data_root <- if(length(args)>0) normalizePath(args[1]) else "../data/raw_dam"
out <- if(length(args)>1) args[2] else "profile_phenotype_discovery"
batch_regex <- if(length(args)>2) args[3] else ".*"
freeze <- "publication_analysis/minute_rebuild"
dir.create(file.path(out,"cache"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,"data"),showWarnings=FALSE)
qc_path <- file.path(freeze,"fly_qc_all_batches.csv")
qc <- fread(qc_path)
stopifnot(!anyDuplicated(qc$id),all(qc$Sex=="M"),all(qc$Temperature=="21.5C"))
keep <- qc[retained_by_existing_gigem==TRUE]
keep[,file:=vapply(strsplit(id,"|",fixed=TRUE),`[`,character(1),2L)]
win <- fread(file.path(freeze,"window_metrics_all_batches.csv"),select=c("id","day","window","observed_bins","sleep_minutes_observed"))
win <- win[window=="ZT00_24"]
batches <- sort(unique(keep$Batch)); batches <- batches[grepl(batch_regex,batches)]
manifest <- list(); source_rows <- list()
event_defs <- data.table(
  event=c("pre_dawn3","pre_dawn6","pre_dusk3","pre_dusk6","pre_dusk30","post_dusk30","post_dawn30","siesta","dark_core"),
  start=c(1260,1080,300,120,450,480,0,120,720),end=c(1440,1440,480,480,480,510,30,360,1200))
for(bi in seq_along(batches)) {
  batch <- batches[bi]; message(sprintf("[%d/%d] %s",bi,length(batches),batch))
  meta <- copy(keep[Batch==batch]); raw_paths <- file.path(data_root,batch,unique(meta$file))
  stopifnot(all(file.exists(raw_paths)))
  inputs <- c(qc_path,file.path(freeze,"window_metrics_all_batches.csv"),
              "scripts/profile_phenotype_helpers.R","scripts/extract_profile_phenotypes.R",raw_paths)
  sig <- unname(tools::md5sum(inputs))
  source_rows[[batch]] <- data.table(batch=batch,path=normalizePath(inputs),md5=sig)
  cache_path <- file.path(out,"cache",paste0(batch,".rds"))
  cached <- if(file.exists(cache_path)) readRDS(cache_path) else NULL
  if(!is.null(cached) && identical(sig,cached$signature) && identical(profile_version,cached$version)) {
    manifest[[batch]] <- copy(cached$manifest); manifest[[batch]][,cached:=TRUE]
    message("  using checksum-validated cache"); next
  }
  loading <- meta[,.(file,region_id,status,start_datetime,stop_datetime,monitor)]
  loading[,start_datetime:=format(start_datetime-600,"%Y-%m-%d %H:%M:%S",tz="UTC")]
  loading[,stop_datetime:=format(stop_datetime,"%Y-%m-%d %H:%M:%S",tz="UTC")]
  linked <- link_dam_metadata(loading,result_dir=file.path(data_root,batch))
  loaded <- load_dam(linked)
  lm <- as.data.table(behavr::meta(loaded))[,.(loaded_id=as.character(id),monitor=as.character(monitor),region_id)]
  idmap <- merge(lm,meta[,.(id,monitor,region_id)],by=c("monitor","region_id"))
  stopifnot(nrow(idmap)==nrow(meta),!anyDuplicated(idmap$loaded_id))
  lu <- setNames(idmap$id,idmap$loaded_id)
  d <- as.data.table(loaded)[,.(id=unname(lu[as.character(id)]),time=as.numeric(t)-600,activity=as.numeric(activity))]
  stopifnot(!anyNA(d$id),all(is.finite(d$activity)),all(d$activity>=0))
  err <- max(abs(d$time/60-round(d$time/60)))
  # Historical DAM files occasionally contain second-level timestamp jitter.
  # Match the publication's nearest-minute mapping and enforce its observed bound.
  stopifnot(err*60 <= max(meta$max_minute_error_seconds)+1e-5, err < .5)
  d[,minute:=as.integer(round(time/60))]; d[,time:=NULL]
  collisions <- d[,.N,by=.(id,minute)][N>1L]
  raw_duplicate_rows <- sum(collisions$N-1L)
  # The frozen pipeline sums count rows mapping to the same nearest minute.
  # Preserve that rule, expose collisions, and require exact downstream agreement.
  if(nrow(collisions)) d <- d[,.(activity=sum(activity)),by=.(id,minute)]
  setorder(d,id,minute)
  d[,moving:=as.numeric(activity>0)]
  d[,run:=cumsum(is.na(shift(minute)) | minute-shift(minute)!=1L | moving!=shift(moving)),by=id]
  d[,asleep:=as.numeric(moving==0 & .N>=5L),by=.(id,run)]
  d <- d[minute>=0 & minute<8640]
  d[,`:=`(day=minute%/%1440L+1L,zt=minute%%1440L,bin=minute%/%30L)]
  check <- d[,.(raw_bins=.N,raw_sleep=sum(asleep)),by=.(id,day)]
  check <- merge(check,win[id %chin% meta$id],by=c("id","day"))
  stopifnot(nrow(check)==nrow(meta)*6L)
  check[,`:=`(coverage_difference=raw_bins-observed_bins,sleep_difference=raw_sleep-sleep_minutes_observed)]
  # Fail on source drift instead of silently mixing live raw data with the freeze.
  stopifnot(all(check$coverage_difference==0),all(check$sleep_difference==0))
  bins <- d[,.(observed=.N,moving=mean(moving),sleep=mean(asleep),activity_rate=mean(activity)),by=.(id,day,bin)]
  grid <- CJ(id=meta$id,bin=0:287)
  bins <- merge(grid,bins,by=c("id","bin"),all.x=TRUE)
  bins[,`:=`(day=bin%/%48L+1L,slot=bin%%48L,observed=fcoalesce(observed,0L))]
  bins[observed<29L,c("moving","sleep","activity_rate"):=NA_real_]
  bins <- merge(bins,meta[,.(id,Batch,Genotype,monitor,region_id)],by="id")
  setorder(bins,id,day,slot)
  daily <- bins[,as.list(profile_features(moving)),by=.(id,day)]
  events <- rbindlist(lapply(seq_len(nrow(event_defs)),function(j) {
    e <- event_defs[j]
    z <- d[zt>=e$start & zt<e$end,.(observed=.N,moving_minutes=sum(moving),rate=mean(moving)),by=.(id,day)]
    z[,event:=e$event]; z[observed<.95*(e$end-e$start),c("moving_minutes","rate"):=NA_real_]; z
  }))
  ev <- dcast(events,id+day~event,value.var=c("moving_minutes","rate"))
  ev[,DawnAnticipation:=fifelse(moving_minutes_pre_dawn6>=10,moving_minutes_pre_dawn3/moving_minutes_pre_dawn6,NA_real_)]
  ev[,DuskAnticipation:=fifelse(moving_minutes_pre_dusk6>=10,moving_minutes_pre_dusk3/moving_minutes_pre_dusk6,NA_real_)]
  ev[,DuskResponse:=rate_post_dusk30-rate_pre_dusk30]
  # A dawn event belongs to the day starting at that dawn; D1 pre-ZT0 is excluded.
  pre <- d[zt>=1410,.(n=.N,pre=mean(moving)),by=.(id,day)]
  pre[n<29L,pre:=NA_real_]; pre[,day:=day+1L]
  ev <- merge(ev,pre[,.(id,day,pre)],by=c("id","day"),all.x=TRUE)
  ev[,DawnResponse:=rate_post_dawn30-pre]
  ev[,`:=`(SiestaMoving=rate_siesta,DarkCoreMoving=rate_dark_core)]
  daily <- merge(daily,ev[,.(id,day,DawnAnticipation,DuskAnticipation,DuskResponse,DawnResponse,SiestaMoving,DarkCoreMoving)],by=c("id","day"),all=TRUE)
  daily <- merge(daily,meta[,.(id,Batch,Genotype,monitor,region_id)],by="id")
  stability <- bins[day>=2,as.list(stability_features(matrix(moving,nrow=5,byrow=TRUE))),by=id]
  stability <- merge(stability,meta[,.(id,Batch,Genotype,monitor,region_id)],by="id")
  row <- data.table(Batch=batch,planned=qc[Batch==batch,.N],retained=nrow(meta),
                    raw_minutes=nrow(d),invalid_halfhour_bins=sum(bins$observed<29L),
                    coverage_max_error=max(abs(check$coverage_difference)),sleep_max_error=max(abs(check$sleep_difference)),raw_duplicate_rows=raw_duplicate_rows,cached=FALSE)
  saveRDS(list(version=profile_version,signature=sig,manifest=row,bins=bins,daily=daily,stability=stability,checks=check,collisions=collisions),cache_path,compress="gzip")
  manifest[[batch]] <- row
  fwrite(rbindlist(manifest),file.path(out,"extraction_manifest.tsv"),sep="\t")
  message(sprintf("  %d retained flies; exact sleep/coverage reconciliation passed",nrow(meta)))
  rm(d,loaded,bins,events,daily); invisible(gc())
}
fwrite(rbindlist(manifest),file.path(out,"extraction_manifest.tsv"),sep="\t")
fwrite(rbindlist(source_rows),file.path(out,"input_checksums.tsv"),sep="\t")
message("Extraction complete: ",normalizePath(out))
