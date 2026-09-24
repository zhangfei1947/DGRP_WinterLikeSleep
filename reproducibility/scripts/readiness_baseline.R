#!/usr/bin/env Rscript
# Read-only sensitivity analyses against the frozen six-day sleep phenotypes.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1]) else normalizePath(".")
out <- file.path(root,"manuscript_readiness","baseline")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source(file.path(root,"scripts","profile_phenotype_helpers.R"))
set.seed(20260910)
nsplit <- 200L
write_tab <- function(x,name) fwrite(x,file.path(out,name),sep="\t",na="NA")
source_daily <- file.path(root,"profile_phenotype_discovery/data/individual_daily_features.rds")
source_windows <- file.path(root,"publication_analysis/minute_rebuild/window_metrics_all_batches.csv")
d <- as.data.table(readRDS(source_daily))
stopifnot(!anyDuplicated(d,by=c("id","day")),all(d$day %in% 1:6))
ids <- unique(d[,.(id,Batch,Genotype)])
stopifnot(nrow(ids)==uniqueN(d$id))
stopifnot(d[Batch=="Batch13_DGRP" & day==1,sum(is.finite(LightSleep))]==0L)

# Strict day requirements: no imputation and no D2-for-D1 substitution.
metrics <- c("DarkSleep","LightSleep","TotalSleep")
features <- list(); paired <- list()
for(m in metrics) {
  z <- dcast(d,id+Batch+Genotype~day,value.var=m)
  M <- as.matrix(z[,as.character(1:6),with=FALSE])
  v <- list(D6_D1=M[,6]-M[,1],Late2_Early2=rowMeans(M[,5:6])-rowMeans(M[,1:2]),
            D6_D2=M[,6]-M[,2],Late2_D2D3=rowMeans(M[,5:6])-rowMeans(M[,2:3]))
  for(k in names(v)) features[[paste(m,k)]] <- cbind(z[,.(id,Batch,Genotype)],data.table(metric=m,contrast=k,value=v[[k]],anchor_metric=m))
  paired[[m]] <- cbind(z[,.(id,Batch,Genotype)],data.table(metric=m,D1=M[,1],D6=M[,6],change=M[,6]-M[,1]))
}
w <- fread(source_windows,select=c("id","Batch","Genotype","day","window","metric_valid","sleep_minutes_scaled"))
w[metric_valid==FALSE,sleep_minutes_scaled:=NA_real_]
segments <- c("ZT00_04","ZT04_08","ZT08_12","ZT12_16","ZT16_20","ZT20_24")
for(s in segments) {
  z <- dcast(w[window==s],id+Batch+Genotype~day,value.var="sleep_minutes_scaled")
  features[[s]] <- z[,.(id,Batch,Genotype,metric=s,contrast="D6_D1",value=`6`-`1`,anchor_metric=if(s %chin% segments[1:2])"LightSleep" else "DarkSleep")]
}
individual <- rbindlist(features)
saveRDS(individual,file.path(out,"individual_sleep_sensitivities.rds"))
lb <- individual[,.(value=finite_median(value),n_flies=sum(is.finite(value))),by=.(Batch,Genotype,metric,contrast,anchor_metric)]
write_tab(lb,"line_batch_sleep_sensitivities.tsv")
line_rows <- list(); quality <- list(); pairs_out <- list(); coupling <- list(); identity <- list(); split_out <- list(); split_design <- list()
for(k in unique(paste(lb$metric,lb$contrast))) {
  z <- lb[paste(metric,contrast)==k]
  fit <- batch_adjust(z)
  eligible <- z[Genotype!="CS" & n_flies>=8 & is.finite(value),.(value=mean(value),n_batches=.N,n_flies=sum(n_flies)),by=Genotype]
  raw <- copy(eligible)[,estimator:="raw_equal_batch_median"]
  adj <- if(is.null(fit)) copy(raw)[0] else merge(fit$lines,eligible[,.(Genotype,n_batches,n_flies)],by="Genotype")[,estimator:="batchFE_median"]
  line_rows[[k]] <- rbind(raw,adj)[,`:=`(metric=z$metric[1],contrast=z$contrast[1],anchor_metric=z$anchor_metric[1])]
  p <- replicate_pairs(z)
  gi <- split(seq_len(nrow(p)),p$Genotype)
  boot <- if(length(gi)>=8) replicate(200,{ix<-unlist(gi[sample(names(gi),length(gi),replace=TRUE)]);pair_rho(p[ix])}) else rep(NA_real_,200)
  qq <- if(sum(is.finite(boot))>20) quantile(boot,c(.025,.975),na.rm=TRUE) else c(NA_real_,NA_real_)
  quality[[k]] <- data.table(metric=z$metric[1],contrast=z$contrast[1],n_lines=nrow(eligible),n_line_batches=nrow(z[Genotype!="CS" & n_flies>=8 & is.finite(value)]),
    n_flies=sum(z[Genotype!="CS" & n_flies>=8 & is.finite(value),n_flies]),n_repeated_lines=length(gi),n_pairs=nrow(p),repeat_rho=pair_rho(p),
    cluster_bootstrap_low=qq[1],cluster_bootstrap_high=qq[2],batchFE_estimable=!is.null(fit))
  if(nrow(p)) pairs_out[[k]] <- p[,`:=`(metric=z$metric[1],contrast=z$contrast[1])]
}
lines <- rbindlist(line_rows); q <- rbindlist(quality)
write_tab(lines,"line_sleep_sensitivities.tsv");write_tab(q,"time_window_repeatability.tsv")
write_tab(rbindlist(pairs_out),"repeatability_pairs.tsv")
agreement <- list()
for(k in unique(paste(lines$metric,lines$contrast,lines$estimator))) {
  z <- lines[paste(metric,contrast,estimator)==k]; a <- lines[metric==z$anchor_metric[1] & contrast=="D6_D1" & estimator==z$estimator[1]]
  x <- merge(z[,.(Genotype,alternative=value)],a[,.(Genotype,anchor=value)],by="Genotype")
  agreement[[k]] <- data.table(metric=z$metric[1],contrast=z$contrast[1],estimator=z$estimator[1],anchor_metric=z$anchor_metric[1],n_common_lines=nrow(x),
    spearman=safe_cor(x$alternative,x$anchor,"spearman"),pearson=safe_cor(x$alternative,x$anchor),same_sign_fraction=mean(sign(x$alternative)==sign(x$anchor)),
    median_alternative_minutes=median(x$alternative),median_anchor_minutes=median(x$anchor))
}
write_tab(rbindlist(agreement),"time_window_agreement.tsv")
# Match flies, not only lines, so removing the D1 requirement is not confounded
# with adding Batch13. Each pair is reaggregated/recalibrated on its common cohort.
matched <- list()
for(m in metrics) for(ct in c("Late2_Early2","D6_D2","Late2_D2D3")) {
  base <- individual[metric==m & contrast=="D6_D1",.(id,Batch,Genotype,anchor=value)]
  alt <- individual[metric==m & contrast==ct,.(id,alternative=value)]
  zz <- merge(base,alt,by="id")[is.finite(anchor) & is.finite(alternative)]
  cell <- zz[,.(n_flies=.N,anchor=median(anchor),alternative=median(alternative)),by=.(Batch,Genotype)]
  eligible <- cell[Genotype!="CS" & n_flies>=8,unique(Genotype)]
  fit_a <- batch_adjust(cell,response="anchor");fit_b <- batch_adjust(cell,response="alternative")
  aa <- merge(fit_a$lines[Genotype %chin% eligible],fit_b$lines[Genotype %chin% eligible],by="Genotype",suffixes=c("_anchor","_alternative"))
  matched[[paste(m,ct)]] <- data.table(metric=m,contrast=ct,n_common_lines=nrow(aa),n_DGRP_flies=cell[Genotype!="CS" & n_flies>=8,sum(n_flies)],
    n_DGRP_cells=cell[Genotype!="CS" & n_flies>=8,.N],spearman=safe_cor(aa$value_anchor,aa$value_alternative,"spearman"),
    same_sign_fraction=mean(sign(aa$value_anchor)==sign(aa$value_alternative)))
}
write_tab(rbindlist(matched),"time_window_common_fly_cohort_agreement.tsv")

# Full-sample associations use paired fly subsets so missingness cannot alter denominators.
# Mean-only identity is kept separate: medians of differences do NOT obey it.
for(m in metrics) {
  z <- paired[[m]][is.finite(D1) & is.finite(D6)]
  cells <- z[,.(n_flies=.N,D1=median(D1),D6=median(D6),change=median(change),mean_D1=mean(D1),mean_D6=mean(D6),mean_change=mean(change)),by=.(Batch,Genotype)]
  eligible <- cells[Genotype!="CS" & n_flies>=8,unique(Genotype)]
  full <- list()
  for(v in c("D1","D6","change")) {
    zr <- cells[,.(Batch,Genotype,n_flies,value=get(v))]
    f <- batch_adjust(zr)
    full[[paste(v,"raw")]] <- zr[Genotype %chin% eligible & n_flies>=8,.(value=mean(value)),by=Genotype][,`:=`(variable=v,estimator="raw_equal_batch_median")]
    if(!is.null(f)) full[[paste(v,"adj")]] <- f$lines[Genotype %chin% eligible][,`:=`(variable=v,estimator="batchFE_median")]
  }
  full <- dcast(rbindlist(full),Genotype+estimator~variable,value.var="value")
  for(e in unique(full$estimator)) {
    x <- full[estimator==e]
    coupling[[paste(m,e)]] <- data.table(metric=m,estimator=e,n_lines=nrow(x),baseline_change_spearman=safe_cor(x$D1,x$change,"spearman"),
      baseline_change_pearson=safe_cor(x$D1,x$change),baseline_D6_spearman=safe_cor(x$D1,x$D6,"spearman"),D1_sd=sd(x$D1),D6_sd=sd(x$D6),change_sd=sd(x$change))
  }
  mn <- cells[Genotype!="CS" & n_flies>=8,.(D1=mean(mean_D1),D6=mean(mean_D6),change=mean(mean_change)),by=Genotype]
  lhs <- cov(mn$D1,mn$change); rhs <- cov(mn$D1,mn$D6)-var(mn$D1)
  stopifnot(abs(lhs-rhs)<1e-7,max(abs(mn$change-(mn$D6-mn$D1)))<1e-7)
  identity[[m]] <- data.table(metric=m,n_lines=nrow(mn),estimator="raw_equal_batch_mean_ONLY",cov_D1_change=lhs,cov_D1_D6_minus_var_D1=rhs,max_mean_identity_error=max(abs(mn$change-(mn$D6-mn$D1))))

  # Independent-fly split diagnostics: retain original >=6 calibration cells,
  # but output lines must have a >=8 full-sample cell (>=4 flies in each half).
  cell_ids <- cells[n_flies>=6,.(Batch,Genotype,n_flies)]
  zz <- merge(z,cell_ids[,.(Batch,Genotype)],by=c("Batch","Genotype"))
  setorder(cell_ids,Batch,Genotype);cell_ids[,cell_id:=.I]
  zz <- merge(zz,cell_ids[,.(Batch,Genotype,cell_id)],by=c("Batch","Genotype")); setorder(zz,cell_id,id)
  ix <- split(seq_len(nrow(zz)),zz$cell_id)
  X <- model.matrix(~factor(Batch)+factor(Genotype),cell_ids)
  Q <- qr(X);stopifnot(Q$rank==ncol(X))
  bc <- grepl("^factor\\(Batch\\)",colnames(X)); be <- X[,bc,drop=FALSE]
  # Calibrated predictions at observed cells minus batch contributions, then equal-batch centre.
  eligible_cell <- cell_ids$Genotype!="CS" & cell_ids$n_flies>=8
  split_design[[m]] <- cell_ids[,.(metric=m,Batch,Genotype,n_full=n_flies,n_half_A=n_flies%/%2L,n_half_B=n_flies-n_flies%/%2L,calibration_cell=TRUE,output_cell=eligible_cell)]
  aggregate_lines <- function(Y,adjust=FALSE) {
    if(adjust) {
      co <- qr.coef(Q,Y);batch_eff <- be %*% co[bc,,drop=FALSE]
      Y <- X %*% co - batch_eff
      # Global intercept shift does not affect any correlation.
    }
    yy <- data.table(Genotype=cell_ids$Genotype[eligible_cell],Y[eligible_cell,,drop=FALSE])
    yy[,lapply(.SD,mean),by=Genotype]
  }
  # Independent implementation of fast multivariate calibration checked against original helper.
  Y0 <- as.matrix(cells[match(paste(cell_ids$Batch,cell_ids$Genotype),paste(cells$Batch,cells$Genotype)),.(D1,change)])
  chk <- aggregate_lines(Y0,TRUE)
  ref <- batch_adjust(cbind(cell_ids[,.(Batch,Genotype,n_flies)],data.table(value=Y0[,1])))$lines
  ref <- merge(chk[,.(Genotype,value_fast=D1)],ref,by="Genotype")
  stopifnot(sd(ref$value_fast-ref$value)<1e-7)
  for(b in seq_len(nsplit)) {
    Y <- t(vapply(ix,function(ii) {
      ia <- sample(ii,length(ii)%/%2L);ib <- setdiff(ii,ia)
      c(baseline_A=median(zz$D1[ia]),change_A=median(zz$change[ia]),change_B=median(zz$change[ib]))
    },numeric(3)))
    for(adj in c(FALSE,TRUE)) {
      g <- aggregate_lines(Y,adj)
      split_out[[paste(m,b,adj)]] <- data.table(metric=m,iteration=b,estimator=if(adj)"batchFE_median" else "raw_equal_batch_median",n_lines=nrow(g),
        shared_half_spearman=safe_cor(g$baseline_A,g$change_A,"spearman"),independent_half_spearman=safe_cor(g$baseline_A,g$change_B,"spearman"),
        shared_half_pearson=safe_cor(g$baseline_A,g$change_A),independent_half_pearson=safe_cor(g$baseline_A,g$change_B))
    }
  }
}
spl <- rbindlist(split_out)
ss <- melt(spl,id.vars=c("metric","iteration","estimator","n_lines"),variable.name="association",value.name="correlation")
sum_s <- ss[,.(n_splits=.N,n_lines=min(n_lines),median=median(correlation),split_q025=quantile(correlation,.025),split_q975=quantile(correlation,.975)),by=.(metric,estimator,association)]
write_tab(rbindlist(coupling),"baseline_change_associations.tsv")
write_tab(rbindlist(identity),"mean_coupling_identity_check.tsv")
write_tab(spl,"independent_fly_split_iterations.tsv");write_tab(sum_s,"independent_fly_split_summary.tsv")
write_tab(spl[,.(n_splits=.N,median_independent_minus_shared=median(independent_half_spearman-shared_half_spearman),
  split_q025=quantile(independent_half_spearman-shared_half_spearman,.025),split_q975=quantile(independent_half_spearman-shared_half_spearman,.975)),by=.(metric,estimator)],
  "independent_fly_split_paired_difference.tsv")
write_tab(rbindlist(split_design),"independent_fly_split_design.tsv")

# Reconcile original individual endpoints rather than merely reusing a downstream table.
original <- as.data.table(readRDS(file.path(root,"profile_phenotype_discovery/data/individual_phenotypes_long.rds")))
v <- individual[metric %chin% metrics & contrast %chin% c("D6_D1","Late2_Early2")]
v[,phenotype:=paste(metric,contrast,sep="__")]
v <- merge(v[,.(id,phenotype,recomputed=value)],original[,.(id,phenotype,original=value)],by=c("id","phenotype"))
validation <- v[,.(n=.N,missingness_disagreements=sum(is.finite(recomputed)!=is.finite(original)),max_abs_difference=max(abs(recomputed-original),na.rm=TRUE)),by=phenotype]
stopifnot(all(validation$missingness_disagreements==0),all(validation$max_abs_difference<1e-8))
write_tab(validation,"original_endpoint_reconciliation.tsv")
dd <- melt(d,id.vars=c("id","day"),measure.vars=metrics,variable.name="metric",value.name="daily_value")
ww <- w[window %chin% c("ZT00_08","ZT08_24","ZT00_24")]
ww[,metric:=fcase(window=="ZT00_08","LightSleep",window=="ZT08_24","DarkSleep",default="TotalSleep")]
dd <- merge(dd,ww[,.(id,day,metric,frozen_value=sleep_minutes_scaled)],by=c("id","day","metric"))
check_frozen <- dd[,.(n=.N,missingness_disagreements=sum(is.finite(daily_value)!=is.finite(frozen_value)),max_abs_difference=max(abs(daily_value-frozen_value),na.rm=TRUE)),by=metric]
stopifnot(all(check_frozen$n==nrow(d)),all(check_frozen$missingness_disagreements==0),all(check_frozen$max_abs_difference<1e-8))
write_tab(check_frozen,"frozen_daily_reconciliation.tsv")
write_tab(d[,.(retained_flies=uniqueN(id),valid_day_sleep=sum(is.finite(LightSleep)),valid_night_sleep=sum(is.finite(DarkSleep)),valid_total_sleep=sum(is.finite(TotalSleep))),by=.(Batch,day)],"daily_coverage_audit.tsv")
manifest <- c(source_daily,source_windows,file.path(root,"scripts/profile_phenotype_helpers.R"),file.path(root,"scripts/readiness_baseline.R"))
write_tab(data.table(path=manifest,md5=unname(tools::md5sum(manifest))),"input_manifest.tsv")
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
message("Baseline and time-window checks complete: ",out)
print(q[metric %chin% metrics]); print(sum_s[estimator=="batchFE_median" & grepl("spearman",association)])
