#!/usr/bin/env Rscript
# Rscript scripts/analyze_profile_phenotypes.R [output_root]
# Descriptive discovery and QA only; no SNP data are tested.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ggrepel);library(lme4)})
source("scripts/profile_phenotype_helpers.R")
setDTthreads(2); set.seed(20260909)
args <- commandArgs(TRUE); out <- if(length(args)) args[1] else "profile_phenotype_discovery"
freeze <- "publication_analysis/minute_rebuild"
for(s in c("data","qa","figures","pca")) dir.create(file.path(out,s),recursive=TRUE,showWarnings=FALSE)
write_tsv <- function(x,p) fwrite(x,file.path(out,p),sep="\t",na="NA")
blue <- "#3174A1"; orange <- "#C77C36"; ink <- "#30343B"
theme_set(theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),plot.title=element_text(colour=ink),strip.background=element_rect(fill="#F2F2F2")))
save_plot <- function(p,name,w=11,h=7) {
  ggsave(file.path(out,"figures",paste0(name,".png")),p,width=w,height=h,dpi=180,bg="white",limitsize=FALSE)
  # Base PDF does not require optional XQuartz/Cairo libraries on macOS.
  ggsave(file.path(out,"figures",paste0(name,".pdf")),p,width=w,height=h,device=grDevices::pdf,useDingbats=FALSE,limitsize=FALSE)
}
manifest <- fread(file.path(out,"extraction_manifest.tsv"))
qc <- fread(file.path(freeze,"fly_qc_all_batches.csv"))
stopifnot(setequal(manifest$Batch,unique(qc$Batch)),sum(manifest$retained)==sum(qc$retained_by_existing_gigem))
caches <- lapply(manifest$Batch,function(b) readRDS(file.path(out,"cache",paste0(b,".rds"))))
bins <- rbindlist(lapply(caches,`[[`,"bins")); daily_new <- rbindlist(lapply(caches,`[[`,"daily"))
stable <- rbindlist(lapply(caches,`[[`,"stability")); checks <- rbindlist(lapply(caches,`[[`,"checks"))
collisions <- rbindlist(lapply(caches,`[[`,"collisions"))
stopifnot(all(checks$sleep_difference==0),all(checks$coverage_difference==0))
rm(caches); invisible(gc())
write_tsv(checks,"qa/raw_sleep_reconciliation.tsv")
write_tsv(collisions,"qa/raw_minute_collisions.tsv")
saveRDS(bins,file.path(out,"data","individual_profiles_30min.rds"))
ids <- unique(bins[,.(id,Batch,Genotype,monitor,region_id)])
write_tsv(qc[,.(Batch,Genotype,monitor,region_id,id,retained_by_existing_gigem,exclusion_reason)],"qa/cohort_membership.tsv")
retention<-qc[,.(planned=.N,retained=sum(retained_by_existing_gigem)),by=.(Batch,Genotype)]
retention[,retained_fraction:=retained/planned]
write_tsv(retention,"qa/cohort_retention_by_line_batch.tsv")

# Daily scalar metrics, inherited from the publication freeze.
wm <- fread(file.path(freeze,"window_metrics_all_batches.csv"))
parts <- wm[window %chin% c("ZT00_08","ZT08_24","ZT00_24")]
parts[,phase:=fcase(window=="ZT00_08","Light",window=="ZT08_24","Dark",default="Total")]
parts[,activity_log:=fifelse(metric_valid,log1p(total_activity/observed_bins),NA_real_)]
parts[,sleep_fraction:=fifelse(metric_valid,sleep_fraction_observed,NA_real_)]
wide <- dcast(parts,id+day~phase,value.var=c("sleep_minutes_scaled","sleep_fraction","activity_log","p_wake_logit","p_doze_logit","fragmentation_log1p","longest_bout_log1p","waking_activity_log1p"))
daily <- wide[,.(id,day,TotalSleep=sleep_minutes_scaled_Total,LightSleep=sleep_minutes_scaled_Light,DarkSleep=sleep_minutes_scaled_Dark,
    SleepAllocation=sleep_fraction_Dark-sleep_fraction_Light,
    ActivityLight=activity_log_Light,ActivityDark=activity_log_Dark,ActivityContrast=activity_log_Dark-activity_log_Light,
    PWakeDark=p_wake_logit_Dark,PDozeDark=p_doze_logit_Dark,FragmentationDark=fragmentation_log1p_Dark,
    LongestBoutDark=longest_bout_log1p_Dark,WakingActivityDark=waking_activity_log1p_Dark,
    PWakeLight=p_wake_logit_Light,PDozeLight=p_doze_logit_Light,FragmentationLight=fragmentation_log1p_Light)]
move <- bins[,.(MovingLight=if(sum(observed[slot<16])>=456) weighted.mean(moving[slot<16],observed[slot<16],na.rm=TRUE) else NA_real_,
               MovingDark=if(sum(observed[slot>=16])>=912) weighted.mean(moving[slot>=16],observed[slot>=16],na.rm=TRUE) else NA_real_),by=.(id,day)]
daily <- merge(daily,move,by=c("id","day"))
daily <- merge(daily,daily_new,by=c("id","day"))
# Keep dawn peak windows contiguous in experimental time: previous evening plus
# the current morning. Folding one changing day at midnight mixes different events.
dawn <- bins[,{
  M <- matrix(moving,nrow=6,byrow=TRUE)
  rbindlist(lapply(1:6,function(k){
    f<-c(MorningPeakOffset=NA_real_,MorningPeakAmplitude=NA_real_)
    if(k>1L){y<-M[k,];y[42:48]<-M[k-1L,42:48];f<-profile_features(y)[names(f)]}
    data.table(day=k,MorningPeakOffset=unname(f[1]),MorningPeakAmplitude=unname(f[2]))
  }))
},by=id]
daily[,c("MorningPeakOffset","MorningPeakAmplitude"):=NULL]
daily<-merge(daily,dawn,by=c("id","day"))
saveRDS(daily,file.path(out,"data","individual_daily_features.rds"))
meta_cols <- c("id","day","Batch","Genotype","monitor","region_id")
metrics <- setdiff(names(daily),meta_cols)
descriptions <- c(
 TotalSleep="24-hour sleep minutes; window coverage >=95%",LightSleep="ZT0-8 sleep minutes",DarkSleep="ZT8-24 sleep minutes",
 SleepAllocation="Dark sleep fraction minus light sleep fraction; window durations normalized",
 ActivityLight="log1p(beam crossings per observed minute), ZT0-8",ActivityDark="log1p(beam crossings per observed minute), ZT8-24",
 ActivityContrast="log1p(dark activity/min) minus log1p(light activity/min); not a raw total ratio",
 MovingLight="Fraction of observed minutes with >=1 beam crossing, ZT0-8",MovingDark="Fraction of observed minutes with >=1 beam crossing, ZT8-24",
 PWakeDark="Corrected log odds of inactive-to-active transitions, ZT8-24; >=30 opportunities",
 PDozeDark="Corrected log odds of active-to-inactive transitions, ZT8-24; >=30 opportunities",
 PWakeLight="Corrected log odds of inactive-to-active transitions, ZT0-8; >=30 opportunities",
 PDozeLight="Corrected log odds of active-to-inactive transitions, ZT0-8; >=30 opportunities",
 FragmentationDark="log1p(bouts per sleep hour), ZT8-24; >=30 sleep minutes; bouts clipped at window boundary",
 FragmentationLight="log1p(bouts per sleep hour), ZT0-8; >=30 sleep minutes; bouts clipped at window boundary",
 LongestBoutDark="log1p(longest sleep bout minutes), clipped to ZT8-24",WakingActivityDark="log1p(beam crossings per awake minute), ZT8-24; >=30 awake minutes",
 Harmonic24Amplitude="24-hour harmonic amplitude from 48-bin moving profile",Harmonic12Amplitude="12-hour harmonic amplitude from 48-bin moving profile",
 TwoHarmonicR2="Variance explained by fixed 24h + 12h harmonics in daily moving profile; LD descriptive only",
 MorningPeakOffset="Local smoothed moving peak offset from dawn; previous-day ZT21-24 plus current ZT0-3; range >=0.10, unique interior maximum",
 EveningPeakOffset="Local smoothed moving peak offset in hours from dusk ZT8, searched ZT5-11; range >=0.10, unique interior maximum",
 MorningPeakAmplitude="Max-minus-min within contiguous previous-day ZT21-24 and current ZT0-3 around dawn",EveningPeakAmplitude="Max-minus-min within ZT5-11",
 DawnAnticipation="Moving minutes ZT21-24 / ZT18-24; >=10 active denominator minutes; relates to next dawn",
 DuskAnticipation="Moving minutes ZT5-8 / ZT2-8; >=10 active denominator minutes",
 DawnResponse="Moving probability ZT0-0.5 minus previous-day ZT23.5-24; D1 unavailable",
 DuskResponse="Moving probability ZT8-8.5 minus ZT7.5-8",
 SiestaMoving="Moving probability ZT2-6",DarkCoreMoving="Moving probability ZT12-20",
 LD24hLagCorrelation="Pearson correlation of adjacent D2-D6 daily moving profiles after centering each day; 48-bin lag",
 LDInterdailyStability="Five times sum of squared mean-profile deviations divided by total D2-D6 squared deviations",
 LDIntradailyVariability="Mean squared successive 30-min difference / population variance, D2-D6",
 OddEvenProfileCorrelation="Correlation of average D2/D4/D6 and D3/D5 moving profiles")
family_for <- function(x) fcase(grepl("P[WD]|Fragmentation|Bout|Waking",x),"sleep_architecture",
  grepl("Peak|Harmonic|TwoHarmonic|LD|OddEven",x),"LD_profile",
  grepl("Anticipation|Response|Siesta|DarkCore",x),"LD_transition",
  grepl("Sleep",x),"sleep_amount_allocation",default="activity")
unit_for <- function(x) fcase(x %chin% c("TotalSleep","LightSleep","DarkSleep"),"minutes",
  grepl("Offset",x),"hours relative to LD transition",grepl("PWake|PDoze",x),"log odds",
  grepl("Activity|Fragmentation|Bout",x),"log1p scale",default="dimensionless")
dictionary <- list(); indiv <- list()
for(m in metrics) {
  dw <- dcast(daily,id~day,value.var=m)
  M <- as.matrix(dw[,as.character(1:6),with=FALSE]); state <- rowMeans(M[,2:6,drop=FALSE])
  phase <- grepl("PeakOffset",m)
  if(phase) {state <- rowMeans(M[,2:6,drop=FALSE],na.rm=TRUE);state[rowSums(is.finite(M[,2:6,drop=FALSE]))<3] <- NA_real_}
  values <- list(State_D2_D6=state,D6_D1=M[,6]-M[,1],Late2_Early2=rowMeans(M[,5:6,drop=FALSE])-rowMeans(M[,1:2,drop=FALSE]))
  if(m %chin% c("DawnResponse","MorningPeakOffset","MorningPeakAmplitude")) values <- values["State_D2_D6"]
  for(cn in names(values)) {
    pn <- paste(m,cn,sep="__")
    indiv[[pn]] <- data.table(id=dw$id,phenotype=pn,value=values[[cn]])
    dictionary[[pn]] <- data.table(phenotype=pn,metric=m,contrast=cn,family=family_for(m),unit=unit_for(m),definition=unname(descriptions[m]),
      aggregation=if(phase && cn=="State_D2_D6") "Mean of >=3 identifiable D2-D6 daily peak offsets" else "Strict within-fly contrast of daily values; all specified days valid",
      eligibility=if(phase) "conditional on identifiable peak; descriptive / validation required" else "finite values and >=8 flies per line-batch",
      selection_class=if(phase) "conditional_phase" else "continuous")
  }
}
for(m in setdiff(names(stable),c("id","Batch","Genotype","monitor","region_id"))) {
  pn <- paste(m,"State_D2_D6",sep="__")
  indiv[[pn]] <- stable[,.(id,phenotype=pn,value=get(m))]
  dictionary[[pn]] <- data.table(phenotype=pn,metric=m,contrast="State_D2_D6",family="LD_profile",unit="dimensionless",definition=unname(descriptions[m]),aggregation="Computed per fly over all 5 complete days",eligibility="240 valid half-hour bins; nonconstant profile; >=8 flies per line-batch",selection_class="continuous")
}

# Profile PCA: fit one common basis per signal and representation, equal weight per
# genotype; exclude CS from training, retain CS for projection/QA. Shape means
# within-profile mean subtraction only (amplitude remains), no noisy unit-variance scaling.
message("Building profile PCA and population profiles")
pca_scores <- list(); loadings <- list(); pca_variance <- list(); pca_stability <- list(); pca_checks <- list(); representatives <- list()
profile_means <- list(); pca_objects <- list(); training_profiles <- list()
for(signal in c("moving","sleep")) {
  for(rep in c("raw","centered_shape","change")) {
    b <- bins[,.(value=if(rep=="change") {
      a<-get(signal)[day %in% 5:6]; e<-get(signal)[day %in% 1:2]; strict_mean(a,2)-strict_mean(e,2)
    } else strict_mean(get(signal)[day>=2],5)),by=.(id,Batch,Genotype,slot)]
    if(rep=="centered_shape") b[,value:=value-mean(value),by=id]
    bw <- dcast(b,id+Batch+Genotype~slot,value.var="value")
    X <- as.matrix(bw[,as.character(0:47),with=FALSE]); ok <- complete.cases(X)
    bw<-bw[ok];X<-X[ok,,drop=FALSE]
    stopifnot(nrow(X)>100,ncol(X)==48)
    if(rep=="raw") profile_means[[signal]] <- b[,.(value=mean(value,na.rm=TRUE),n_flies=sum(is.finite(value))),by=.(Batch,Genotype,slot)][,signal:=signal]
    lbp <- bw[,lapply(.SD,mean),by=.(Batch,Genotype),.SDcols=as.character(0:47)]
    nn <- bw[,.N,by=.(Batch,Genotype)];setnames(nn,"N","n_flies");lbp<-merge(lbp,nn,by=c("Batch","Genotype"))
    adj <- lapply(as.character(0:47),function(k) batch_adjust(lbp,k))
    stopifnot(all(vapply(adj,function(a)!is.null(a),logical(1))))
    eligible <- lbp[Genotype!="CS" & n_flies>=8,unique(Genotype)]
    G <- Reduce(function(a,b) merge(a,b,by="Genotype"),lapply(seq_along(adj),function(k){a<-copy(adj[[k]]$lines);setnames(a,"value",paste0("bin",k));a}))
    G <- G[Genotype %chin% eligible];setorder(G,Genotype)
    GX <- as.matrix(G[,-1]); model <- prcomp(GX,center=TRUE,scale.=FALSE)
    # Deterministic orientation: largest absolute loading is positive.
    for(k in seq_len(ncol(model$rotation))) {sgn<-sign(model$rotation[which.max(abs(model$rotation[,k])),k]); model$rotation[,k]<-model$rotation[,k]*sgn;model$x[,k]<-model$x[,k]*sgn}
    scores <- sweep(X,2,model$center,"-")%*%model$rotation[,1:3,drop=FALSE]
    stem <- paste(signal,rep,sep="_")
    training_profiles[[stem]] <- copy(G)[,model:=stem]
    pca_objects[[stem]] <- model
    prop <- model$sdev^2/sum(model$sdev^2)
    pca_variance[[stem]] <- data.table(model=stem,PC=seq_along(prop),variance_fraction=prop,n_training_lines=nrow(G),n_projected_flies=nrow(X))
    for(k in 1:3) {
      pn<-paste0("PCA_",stem,"_PC",k)
      pca_scores[[pn]]<-data.table(id=bw$id,phenotype=pn,value=scores[,k])
      dictionary[[pn]]<-data.table(phenotype=pn,metric=pn,contrast=if(rep=="change")"Late2_Early2" else "State_D2_D6",family="profile_PCA",unit="PC score",definition=paste("PC",k,"of",signal,rep,"48-bin profiles; basis trained on equal-weight batch-adjusted DGRP line means"),aggregation="Per-fly projection; mean score per line-batch (linear projection); batch-adjusted line score",eligibility="All bins valid for required days; >=8 flies per line-batch",selection_class="PCA_discovery")
      loadings[[paste(stem,k)]]<-data.table(model=stem,PC=paste0("PC",k),zt=seq(.25,23.75,.5),loading=model$rotation[,k],centre=model$center,score_sd=model$sdev[k],variance_fraction=prop[k])
      oo<-order(model$x[,k]);chosen<-unique(c(oo[1:3],tail(oo,3)))
      representatives[[paste(stem,k)]]<-data.table(model=stem,PC=k,Genotype=G$Genotype[chosen],score=model$x[chosen,k])
    }
    # Basis stability on genotype halves: no labels or SNP associations are used.
    ix<-sample(seq_len(nrow(G)));a<-ix[seq_len(floor(length(ix)/2))];bb<-setdiff(ix,a)
    A<-prcomp(GX[a,,drop=FALSE])$rotation[,1:3,drop=FALSE];B<-prcomp(GX[bb,,drop=FALSE])$rotation[,1:3,drop=FALSE]
    sv<-svd(crossprod(A,B))$d
    pca_stability[[stem]]<-data.table(model=stem,subspace_axis=1:3,cos_principal_angle=sv,n_half1=length(a),n_half2=length(bb))
    pca_checks[[stem]]<-data.table(model=stem,orthogonality_max_error=max(abs(crossprod(model$rotation)-diag(ncol(model$rotation)))),shape_mean_max=if(rep=="centered_shape")max(abs(rowMeans(X))) else NA_real_)
  }
}
saveRDS(pca_objects,file.path(out,"pca","models.rds"))
write_tsv(rbindlist(training_profiles),"pca/training_line_profiles.tsv")
write_tsv(rbindlist(pca_variance),"pca/variance_explained.tsv");write_tsv(rbindlist(loadings),"pca/loadings.tsv")
write_tsv(rbindlist(pca_stability),"pca/split_genotype_basis_stability.tsv");write_tsv(rbindlist(pca_checks),"qa/pca_checks.tsv")
write_tsv(rbindlist(representatives),"pca/extreme_line_examples.tsv")
pop <- rbindlist(profile_means);write_tsv(pop,"data/population_profiles_D2_D6.tsv")
feature_grid <- CJ(id=ids$id,phenotype=names(dictionary))
features <- merge(feature_grid,rbindlist(c(indiv,pca_scores)),by=c("id","phenotype"),all.x=TRUE)
features <- merge(features,ids,by="id",all.x=TRUE)
features[!is.finite(value),value:=NA_real_]
dict <- rbindlist(dictionary,fill=TRUE)
stopifnot(!anyDuplicated(dict$phenotype),!anyDuplicated(features,by=c("id","phenotype")))
saveRDS(features,file.path(out,"data","individual_phenotypes_long.rds"))
write_tsv(dict,"phenotype_dictionary.tsv")
lb <- features[,.(value=if(grepl("^PCA_",phenotype[1])) finite_mean(value) else finite_median(value),
                  n_flies=sum(is.finite(value)),n_retained=.N,mad=if(sum(is.finite(value))>1) mad(value,na.rm=TRUE) else NA_real_,
                  q25=if(any(is.finite(value))) quantile(value,.25,na.rm=TRUE) else NA_real_,q75=if(any(is.finite(value)))quantile(value,.75,na.rm=TRUE) else NA_real_),by=.(Batch,Genotype,phenotype)]
lb<-merge(lb,dict[,.(phenotype,metric,contrast,family,selection_class)],by="phenotype")
write_tsv(lb,"data/line_batch_phenotypes.tsv")
availability <- features[,.(n_flies=.N,valid_flies=sum(is.finite(value)),valid_fraction=mean(is.finite(value))),by=.(phenotype,Batch)]
write_tsv(availability,"qa/phenotype_availability_by_batch.tsv")

message("Evaluating repeatability, batch sensitivity and minimum sample thresholds")
qa <- list(); raw_lines<-list(); adj_lines<-list(); pair_all<-list(); min_n_rows<-list(); batch_effects<-list()
for(j in seq_len(nrow(dict))) {
  pn<-dict$phenotype[j];z<-lb[phenotype==pn];p<-replicate_pairs(z);rho<-pair_rho(p)
  if(j%%15==0) message(j,"/",nrow(dict)," phenotypes assessed")
  gp<-unique(p$Genotype)
  pair_ix<-split(seq_len(nrow(p)),p$Genotype)
  boot<-if(length(gp)>=8) replicate(200,{ix<-unlist(pair_ix[sample(gp,length(gp),replace=TRUE)],use.names=FALSE);safe_cor(c(p$x[ix],p$y[ix]),c(p$y[ix],p$x[ix]),"spearman")}) else rep(NA_real_,200)
  ci<-if(sum(is.finite(boot))>20) quantile(boot,c(.025,.975),na.rm=TRUE) else c(NA_real_,NA_real_)
  loo<-if(length(gp)>=8) vapply(gp,function(g)pair_rho(p[Genotype!=g]),numeric(1)) else NA_real_
  raw<-z[Genotype!="CS" & n_flies>=8 & is.finite(value),.(value=mean(value),n_batches=.N,total_flies=sum(n_flies)),by=Genotype]
  raw[,phenotype:=pn]; raw_lines[[pn]]<-raw
  fit<-batch_adjust(z);adj<-if(!is.null(fit)) merge(fit$lines,raw[,.(Genotype,n_batches,total_flies)],by="Genotype") else data.table(Genotype=character(),value=numeric(),n_batches=integer(),total_flies=integer())
  adj[,phenotype:=pn];adj_lines[[pn]]<-adj
  if(!is.null(fit)) batch_effects[[pn]]<-copy(fit$batches)[,phenotype:=pn]
  cross<-merge(raw[,.(Genotype,raw=value)],adj[,.(Genotype,adjusted=value)],by="Genotype")
  csz<-z[Genotype=="CS" & n_flies>=8 & is.finite(value)]
  ss<-sd(raw$value);csratio<-if(nrow(csz)>=3 && is.finite(ss) && ss>0)sd(csz$value)/ss else NA_real_
  var_ratio<-NA_real_;singular<-NA
  zz<-z[n_flies>=8 & Genotype!="CS" & is.finite(value)]
  if(length(gp)>=10 && nrow(zz)>uniqueN(zz$Genotype)+10 && sd(zz$value)>1e-8) {
    mm<-tryCatch(suppressWarnings(lmer(value~Batch+(1|Genotype),data=zz,REML=TRUE,control=lmerControl(check.conv.singular="ignore"))),error=function(e)NULL)
    if(!is.null(mm)){vc<-as.data.frame(VarCorr(mm));var_ratio<-vc$vcov[vc$grp=="Genotype"]/sum(vc$vcov);singular<-isSingular(mm)}
  }
  v16<-z[Genotype!="CS" & n_flies>=16 & is.finite(value),.(value=mean(value)),by=Genotype]
  xx<-merge(raw[,.(Genotype,v8=value)],v16,by="Genotype")
  m16rho<-safe_cor(xx$v8,xx$value,"spearman")
  for(k in c(8,12,16)) {pk<-replicate_pairs(z,k);min_n_rows[[paste(pn,k)]]<-data.table(phenotype=pn,min_n=k,lines=uniqueN(z[n_flies>=k & Genotype!="CS" & is.finite(value),Genotype]),repeated_lines=uniqueN(pk$Genotype),pairs=nrow(pk),symmetric_rho=pair_rho(pk))}
  qa[[pn]]<-data.table(phenotype=pn,eligible_lines=nrow(raw),eligible_line_batches=nrow(zz),valid_fly_fraction=sum(z$n_flies)/sum(z$n_retained),
      repeated_lines=length(gp),batch_pairs=nrow(p),repeatability_rho=rho,rho_ci_low=ci[1],rho_ci_high=ci[2],
      traditional_pair_rho=if(nrow(p)>=8)safe_cor(p$x,p$y,"spearman") else NA_real_,
      max_LOO_rho_change=if(any(is.finite(loo)))max(abs(loo-rho),na.rm=TRUE) else NA_real_,
      raw_batchFE_rho=safe_cor(cross$raw,cross$adjusted,"spearman"),batchFE_available=!is.null(fit),
      CS_batches=nrow(csz),CS_sd_vs_DGRP_sd=csratio,conditional_line_variance_fraction=var_ratio,mixed_model_singular=singular,
      min8_min16_line_rho=m16rho,min16_lines=nrow(v16),raw_sd=sd(raw$value))
  if(nrow(p)){p[,phenotype:=pn];pair_all[[pn]]<-p}
}
quality<-merge(rbindlist(qa),dict,by="phenotype")
quality[,review_status:=fcase(selection_class=="conditional_phase","phase-selected: validate missingness",
    eligible_lines<100 | repeated_lines<20 | !batchFE_available,"insufficient coverage/calibration",
    is.finite(repeatability_rho) & repeatability_rho>=.5 & rho_ci_low>.2 & raw_batchFE_rho>=.7 & min8_min16_line_rho>=.8,"priority candidate",
    is.finite(repeatability_rho) & repeatability_rho>=.35,"secondary / batch-sensitive",default="exploratory / low repeatability")]
setorder(quality,-repeatability_rho)
write_tsv(quality,"qa/phenotype_quality_rank.tsv");write_tsv(rbindlist(min_n_rows),"qa/minimum_n_sensitivity.tsv")
write_tsv(rbindlist(pair_all),"qa/repeated_line_pairs.tsv");write_tsv(rbindlist(batch_effects),"qa/batch_effects.tsv")
raw<-rbindlist(raw_lines);adj<-rbindlist(adj_lines)
write_tsv(raw,"data/line_phenotypes_raw.tsv");write_tsv(adj,"data/line_phenotypes_batchFE.tsv")
rawwide<-dcast(raw,Genotype~phenotype,value.var="value");adjwide<-dcast(adj,Genotype~phenotype,value.var="value")
write_tsv(rawwide,"data/line_matrix_raw.tsv");write_tsv(adjwide,"data/line_matrix_batchFE.tsv")
valid_cols<-names(adjwide)[-1];valid_cols<-valid_cols[vapply(adjwide[,..valid_cols],function(x)sum(is.finite(x))>=50 && sd(x,na.rm=TRUE)>1e-8,logical(1))]
cm<-cor(as.matrix(adjwide[,..valid_cols]),use="pairwise.complete.obs",method="spearman")
cn<-crossprod(is.finite(as.matrix(adjwide[,..valid_cols])))
write_tsv(data.table(phenotype=rownames(cm),as.data.table(cm)),"qa/phenotype_spearman_matrix.tsv")
write_tsv(data.table(phenotype=rownames(cn),as.data.table(cn)),"qa/correlation_pairwise_n.tsv")
redundant<-as.data.table(as.table(cm));setnames(redundant,c("phenotype1","phenotype2","rho"))
redundant<-redundant[as.character(phenotype1)<as.character(phenotype2) & abs(rho)>=.85]
write_tsv(redundant,"qa/redundant_phenotype_pairs_absrho_ge_085.tsv")
pca_cor<-as.data.table(as.table(cm));setnames(pca_cor,c("PCA_phenotype","comparison_phenotype","rho"))
pca_cor<-pca_cor[grepl("^PCA_",PCA_phenotype) & !grepl("^PCA_",comparison_phenotype)]
pca_cor[,paired_lines:=mapply(function(a,b)cn[a,b],PCA_phenotype,comparison_phenotype)]
setorder(pca_cor,PCA_phenotype,-rho)
write_tsv(pca_cor,"pca/scalar_interpretation_correlations.tsv")
# Greedy discovery panel, limited per family and retaining the original endpoint as anchor.
pool<-quality[review_status=="priority candidate" & is.finite(max_LOO_rho_change) & max_LOO_rho_change<.15]
selected<-"DarkSleep__D6_D1";family_counts<-setNames(rep(0,6),c("sleep_architecture","LD_profile","LD_transition","sleep_amount_allocation","profile_PCA","activity"))
family_counts["sleep_amount_allocation"]<-1L
for(pn in pool$phenotype) {
  if(pn %chin% selected || length(selected)>=10 || !pn %in% colnames(cm)) next
  f<-quality[phenotype==pn,family]
  if(f %in% names(family_counts) && family_counts[f]>=2) next
  prev<-intersect(selected,colnames(cm));rho<-cm[pn,prev]
  if(length(rho) && any(abs(rho)>.85,na.rm=TRUE)) next
  selected<-c(selected,pn);if(f%in%names(family_counts))family_counts[f]<-family_counts[f]+1
}
short<-quality[match(selected,phenotype)]
short[,panel_role:=ifelse(phenotype=="DarkSleep__D6_D1","existing endpoint anchor","pilot candidate; validate independently")]
write_tsv(short,"candidate_panel.tsv")

message("Rendering standalone QA figures")
# Chart contract: static ggplot PNG/PDF; two palette roots, direct labels; rank
# plot uses genotype-cluster CIs, heatmap uses paired line counts in source TSV,
# profile plots use 48 ZT bins, 8L:16D annotation, and explicit sample definitions.
show_names<-unique(c(head(quality[selection_class!="conditional_phase" & repeated_lines>=20,phenotype],28),selected))
show<-quality[phenotype %chin% show_names]
show[,label:=gsub("__"," / ",phenotype,fixed=TRUE)];show[,label:=factor(label,levels=rev(label))]
p<-ggplot(show,aes(repeatability_rho,label))+geom_vline(xintercept=c(0,.5),colour="grey75",linetype="dashed")+
  geom_segment(aes(x=rho_ci_low,xend=rho_ci_high,yend=label),colour="grey55")+geom_point(colour=blue,size=2)+
  labs(title="Cross-batch phenotype repeatability",subtitle="Symmetric Spearman correlation; 95% genotype-cluster bootstrap CI; >=8 flies per cell",x="Repeatability (rho)",y=NULL)
save_plot(p,"01_repeatability",12,max(10,.30*nrow(show)+1.5))
sel<-intersect(selected,colnames(cm));cm2<-cm[sel,sel,drop=FALSE];cc<-as.data.table(as.table(cm2));setnames(cc,c("a","b","rho"))
p<-ggplot(cc,aes(a,b,fill=rho))+geom_tile(colour="white")+geom_text(aes(label=sprintf("%.2f",rho)),size=3)+
  scale_fill_gradient2(low=orange,mid="white",high=blue,limits=c(-1,1))+coord_equal()+
  labs(title="Candidate phenotype correlations",subtitle="Spearman correlation of batch-adjusted DGRP line estimates; pairwise complete observations",x=NULL,y=NULL,fill="rho")+
  theme(axis.text.x=element_text(angle=50,hjust=1,size=8),axis.text.y=element_text(size=8))
save_plot(p,"02_candidate_correlations",12,10)
ld<-rbindlist(loadings);ld[,PC:=factor(PC,levels=c("PC1","PC2","PC3"))]
p<-ggplot(ld,aes(zt,loading))+annotate("rect",xmin=8,xmax=24,ymin=-Inf,ymax=Inf,fill="grey93")+geom_hline(yintercept=0,colour="grey65")+
  geom_line(colour=blue,linewidth=.6)+facet_grid(model~PC,scales="free_y")+scale_x_continuous(breaks=c(0,8,16,24),limits=c(0,24))+
  labs(title="Profile PCA loading curves",subtitle="Common bases from equal-weight batch-adjusted line means; shaded ZT8-24 is dark",x="Zeitgeber time (hours)",y="Loading")
save_plot(p,"03_PCA_loadings",12,13)
recon<-rbindlist(lapply(split(ld,interaction(ld$model,ld$PC)),function(z) {
 rbind(z[,.(model,PC,zt,curve="minus 1 SD",value=centre-loading*score_sd)],z[,.(model,PC,zt,curve="plus 1 SD",value=centre+loading*score_sd)])
}))
p<-ggplot(recon,aes(zt,value,colour=curve,linetype=curve))+annotate("rect",xmin=8,xmax=24,ymin=-Inf,ymax=Inf,fill="grey93")+geom_line(linewidth=.6)+
 facet_grid(model~PC,scales="free_y")+scale_colour_manual(values=c(orange,blue))+scale_linetype_manual(values=c("dashed","solid"))+
 scale_x_continuous(breaks=c(0,8,16,24))+labs(title="PCA profile reconstruction",subtitle="Centre +/- one component SD; model illustration, not observed individual probabilities",x="Zeitgeber time (hours)",y="Profile value",colour=NULL,linetype=NULL)+theme(legend.position="top")
save_plot(p,"04_PCA_reconstruction",12,13)
cs<-lb[Genotype=="CS" & phenotype %chin% selected & n_flies>=8 & is.finite(value)]
reference<-raw[,.(med=median(value),spread=sd(value)),by=phenotype]
cs<-merge(cs,reference,by="phenotype");cs[,relative:=(value-med)/spread]
cs[,batch_number:=as.numeric(sub("[a-z]?_DGRP$","",sub("Batch","",Batch)))]
p<-ggplot(cs,aes(batch_number,relative))+geom_hline(yintercept=0,colour="grey65")+geom_point(colour=blue)+geom_line(colour=blue)+
  facet_wrap(~phenotype,ncol=2)+labs(title="CS profiles across experimental batches",subtitle="Raw CS cell estimates relative to the DGRP line median, scaled by DGRP line SD; gaps are missing CS batches",x="Batch number",y="Difference / DGRP SD")
save_plot(p,"05_CS_stability",12,max(7,2.2*ceiling(length(unique(cs$phenotype))/2)))
chosen<-head(short$phenotype,6);pp<-rbindlist(pair_all)[phenotype %chin% chosen]
p<-ggplot(pp,aes(x,y))+geom_abline(slope=1,intercept=0,colour="grey60",linetype="dashed")+geom_point(colour=blue,alpha=.6,size=1.5)+
  geom_text_repel(data=pp[Genotype=="DGRP509"],aes(label=Genotype),max.overlaps=Inf,size=3)+
  facet_wrap(~phenotype,scales="free",ncol=2)+labs(title="Independent-batch phenotype pairs",subtitle="Each point is one pair of batches for the same DGRP line; >=8 flies per line-batch",x="Earlier-sorted batch value",y="Later-sorted batch value")
save_plot(p,"06_replicate_pairs",12,10)
peak<-daily[,.(Morning=mean(is.finite(MorningPeakOffset)),Evening=mean(is.finite(EveningPeakOffset))),by=.(Batch,Genotype)]
write_tsv(peak,"qa/peak_identifiability.tsv")
allcs<-pop[Genotype=="CS" & n_flies>=8]
p<-ggplot(allcs,aes(slot*.5+.25,value,group=Batch))+annotate("rect",xmin=8,xmax=24,ymin=-Inf,ymax=Inf,fill="grey93")+
 geom_line(colour=blue,alpha=.28)+facet_wrap(~signal,ncol=1)+scale_x_continuous(breaks=c(0,8,16,24),limits=c(0,24))+
 labs(title="CS mean D2-D6 profiles",subtitle="One curve per batch; same retained flies throughout D2-D6; shaded ZT8-24 is dark",x="Zeitgeber time (hours)",y="Probability")
save_plot(p,"07_CS_profile_overlays",10,7)
ext<-rbindlist(representatives)[model=="moving_centered_shape" & PC %in% 1:2,unique(Genotype)]
ep<-pop[signal=="moving" & Genotype %chin% ext & n_flies>=8]
ep_n<-unique(ep[,.(Batch,Genotype,n_flies)])[,.(label=paste0(Genotype[1]," (n=",sum(n_flies),", batches=",.N,")")),by=Genotype]
ep<-merge(ep,ep_n,by="Genotype")
p<-ggplot(ep,aes(slot*.5+.25,value,group=Batch))+annotate("rect",xmin=8,xmax=24,ymin=-Inf,ymax=Inf,fill="grey93")+
 geom_line(colour=blue,alpha=.65)+facet_wrap(~label,ncol=3)+scale_x_continuous(breaks=c(0,8,16,24))+
 labs(title="Observed profiles at extremes of moving-profile shape PCs",subtitle="Raw D2-D6 population means; one curve per batch; selected using batch-adjusted PC1/PC2",x="Zeitgeber time (hours)",y="Moving probability")
save_plot(p,"08_shape_PC_extreme_profiles",12,3*ceiling(length(unique(ep$Genotype))/3))

# Compact run summary and machine-readable verification.
legacy<-fread(file.path(freeze,"line_batch_phenotypes_long.csv"))[phenotype=="NightSleep_D6_D1",.(Batch,Genotype,legacy_value=value,legacy_n=n_flies)]
legacy<-merge(legacy,lb[phenotype=="DarkSleep__D6_D1",.(Batch,Genotype,value,n_flies)],by=c("Batch","Genotype"),all=TRUE)
write_tsv(legacy,"qa/legacy_nightsleep_reconciliation.tsv")
validation<-data.table(check=c("frozen_cohort_count","raw_sleep_reconciliation","raw_coverage_reconciliation","feature_key_unique","no_infinite_features","PCA_orthogonality","centered_shape_zero_mean","Batch13_D1_total_missing","legacy_nightsleep_reconciliation"),
  passed=c(nrow(ids)==sum(qc$retained_by_existing_gigem),all(checks$sleep_difference==0),all(checks$coverage_difference==0),!anyDuplicated(features,by=c("id","phenotype")),!any(is.infinite(features$value)),
      all(rbindlist(pca_checks)$orthogonality_max_error<1e-8),all(rbindlist(pca_checks)[is.finite(shape_mean_max),shape_mean_max]<1e-8),
      all(!is.finite(daily[Batch=="Batch13_DGRP" & day==1,TotalSleep])),
      all(legacy$legacy_n==legacy$n_flies) && all(abs(legacy$legacy_value-legacy$value)<1e-8)))
write_tsv(validation,"qa/validation_checks.tsv");stopifnot(all(validation$passed))
counts<-data.table(metric=c("batches","retained_flies","DGRP_lines_before_min_n","phenotypes","priority_candidates_before_redundancy","panel_size"),value=c(nrow(manifest),nrow(ids),uniqueN(ids[Genotype!="CS",Genotype]),nrow(dict),quality[review_status=="priority candidate",.N],nrow(short)))
write_tsv(counts,"run_summary.tsv")
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
code_files<-c("scripts/profile_phenotype_helpers.R","scripts/extract_profile_phenotypes.R","scripts/analyze_profile_phenotypes.R")
write_tsv(data.table(file=normalizePath(code_files),md5=unname(tools::md5sum(code_files))),"analysis_code_checksums.tsv")
fig_files<-list.files(file.path(out,"figures"),pattern="\\.(png|pdf)$",full.names=TRUE)
stopifnot(length(fig_files)==16L,all(file.info(fig_files)$size>1000))
write_tsv(data.table(file=basename(fig_files),bytes=file.info(fig_files)$size,md5=unname(tools::md5sum(fig_files))),"qa/figure_manifest.tsv")
message("Completed profile discovery in ",normalizePath(out))
