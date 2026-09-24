#!/usr/bin/env Rscript
# Rscript scripts/readiness_pca.R [output_dir] [line_B=500] [fly_B=200] [selection_B=500]
# Manuscript sensitivity only: frozen primary phenotype definitions and GWAS unchanged.
suppressPackageStartupMessages(library(data.table))
source("scripts/profile_phenotype_helpers.R")
setDTthreads(2)
args <- commandArgs(TRUE)
out <- if(length(args)) args[1] else "manuscript_readiness/pca"
line_B <- if(length(args)>1) as.integer(args[2]) else 500L
fly_B <- if(length(args)>2) as.integer(args[3]) else 200L
selection_B <- if(length(args)>3) as.integer(args[4]) else 500L
stopifnot(line_B>=200L, fly_B>=100L, selection_B>=200L)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
write_tsv <- function(x,name) fwrite(x,file.path(out,name),sep="\t",na="NA")
src <- "profile_phenotype_discovery"
input_files <- c("scripts/readiness_pca.R", "scripts/profile_phenotype_helpers.R",
  file.path(src,c("data/individual_profiles_30min.rds","pca/models.rds",
  "pca/training_line_profiles.tsv","qa/phenotype_quality_rank.tsv",
  "qa/repeated_line_pairs.tsv","qa/phenotype_spearman_matrix.tsv","candidate_panel.tsv")))
manifest <- data.table(path=input_files, md5=unname(tools::md5sum(input_files)))
stopifnot(!anyNA(manifest$md5))
write_tsv(manifest,"input_manifest.tsv")
writeLines(c("Seeds: line bootstrap 202609101; whole-fly bootstrap 202609102; selection 202609103.",
  sprintf("Replicates: line=%d; within-cell whole-fly=%d; conditional selection=%d.",line_B,fly_B,selection_B),
  "No inferential acceptance gates or primary-phenotype replacements are introduced.",
  "Each PCA fit uses equal DGRP-line weights; CS calibrates batch effects but never trains PCA.",
  "PCA signs align with frozen loadings. Subspace cosines ignore arbitrary signs/rotations.",
  "Line bootstrap holds batch-adjusted profiles fixed; fly bootstrap refits batch correction and PCA with the original design/cohort fixed.",
  "Selection holds original QA pool, correlations, family cap and LOO gate fixed, and resamples genotype clusters jointly across traits for ranking only."),
  file.path(out,"PROTOCOL.txt"))

models <- readRDS(file.path(src,"pca/models.rds"))
training <- fread(file.path(src,"pca/training_line_profiles.tsv"))
model_names <- c("moving_raw","moving_centered_shape")
GX <- lapply(setNames(model_names,model_names), function(stem) {
  z <- training[model==stem]; setorder(z,Genotype)
  m <- as.matrix(z[,paste0("bin",1:48),with=FALSE]); rownames(m)<-z$Genotype; m
})
stopifnot(identical(rownames(GX[[1]]),rownames(GX[[2]])))
metrics <- function(fit, stem, perturbed=NULL, held_out=NULL) {
  ref <- models[[stem]]; g <- GX[[stem]]
  cc <- crossprod(ref$rotation[,1:3], fit$rotation[,1:3])
  sgn <- if(cc[1,1]<0) -1 else 1
  v <- fit$rotation[,1]*sgn
  old_score <- as.vector(sweep(g,2,ref$center,"-") %*% ref$rotation[,1])
  projected <- as.vector(sweep(g,2,fit$center,"-") %*% v)
  sc <- if(is.null(perturbed)) NA_real_ else safe_cor(old_score,
    as.vector(sweep(perturbed,2,fit$center,"-") %*% v),"spearman")
  variance <- fit$sdev^2/sum(fit$sdev^2)
  data.table(PC1_loading_cosine=abs(cc[1,1]),
    frozen_PC2_to_new_PC1_abs_cosine=abs(cc[2,1]),
    PC1_best_matching_frozen_axis=which.max(abs(cc[,1])),
    PC1_fixed_profile_score_rho=safe_cor(old_score,projected,"spearman"),
    PC1_perturbed_profile_score_rho=sc,
    PC1_out_of_bag_score_rho=if(length(held_out)>=8) safe_cor(old_score[held_out],projected[held_out],"spearman") else NA_real_,
    n_out_of_bag=length(held_out),
    first2_subspace_min_cosine=min(svd(cc[1:2,1:2])$d),
    first3_subspace_min_cosine=min(svd(cc)$d),
    variance_fraction_PC1=variance[1],variance_fraction_PC2=variance[2],
    variance_fraction_PC3=variance[3],variance_gap_PC1_PC2=variance[1]-variance[2])
}

message("PCA line-bootstrap refits (frozen batch adjustment)")
set.seed(202609101)
draws <- replicate(line_B,sample.int(nrow(GX[[1]]),replace=TRUE),simplify=FALSE)
results <- list(); loadings <- list(); ri <- 0L; li <- 0L
for(stem in model_names) for(b in seq_len(line_B)) {
  fit <- prcomp(GX[[stem]][draws[[b]],,drop=FALSE],center=TRUE,scale.=FALSE)
  z <- metrics(fit,stem,held_out=setdiff(seq_len(nrow(GX[[stem]])),draws[[b]]))
  z[,`:=`(model=stem,resampling="DGRP_line_bootstrap_frozen_batch",replicate=b)]
  results[[ri<-ri+1L]] <- z
  sgn <- sign(sum(models[[stem]]$rotation[,1]*fit$rotation[,1]))
  loadings[[li<-li+1L]] <- data.table(model=stem,resampling=z$resampling,replicate=b,
    bin=1:48,loading=fit$rotation[,1]*sgn)
}

message("Reconstructing individual D2-D6 profiles and independent multivariate batch calibration")
bins <- readRDS(file.path(src,"data/individual_profiles_30min.rds"))
curves <- bins[day>=2,.(value=if(.N==5L && all(is.finite(moving))) mean(moving) else NA_real_),
  by=.(id,Batch,Genotype,slot)]
wide <- dcast(curves,id+Batch+Genotype~slot,value.var="value")
X <- as.matrix(wide[,as.character(0:47),with=FALSE])
ok <- complete.cases(X); wide<-wide[ok]; X<-X[ok,,drop=FALSE]
rm(bins,curves); invisible(gc())
cells <- unique(wide[,.(Batch,Genotype)]); setorder(cells,Batch,Genotype)
cells[,cell:=.I]
wide <- merge(wide,cells,by=c("Batch","Genotype"),sort=FALSE)
# merge order is not assumed: rebuild X after attaching cell identity.
X <- as.matrix(wide[,as.character(0:47),with=FALSE])
cell_rows <- split(seq_len(nrow(X)),factor(wide$cell,levels=cells$cell))
cells[,n_flies:=lengths(cell_rows)]
fit_cells <- which(cells$n_flies>=6L)
design_data <- copy(cells[fit_cells])
design_data[,`:=`(Batch=factor(Batch),Genotype=factor(Genotype))]
D <- model.matrix(~Batch+Genotype,design_data)
qd <- qr(D); stopifnot(qd$rank==ncol(D))
g_names <- rownames(GX[[1]])
stopifnot(setequal(g_names,cells[Genotype!="CS" & n_flies>=8,unique(Genotype)]))
# EMM contrast averaged over every fitted batch; constructed independently of batch_adjust().
newdata <- CJ(Genotype=g_names,Batch=levels(design_data$Batch))
newdata[,`:=`(Batch=factor(Batch,levels=levels(design_data$Batch)),
  Genotype=factor(Genotype,levels=levels(design_data$Genotype)))]
P <- model.matrix(~Batch+Genotype,newdata)
E <- rowsum(P,group=newdata$Genotype,reorder=TRUE)/nlevels(design_data$Batch)
E <- E[g_names,,drop=FALSE]
map <- E %*% qr.coef(qd,diag(nrow(D)))
cell_means <- t(vapply(cell_rows,function(ix)colMeans(X[ix,,drop=FALSE]),numeric(48)))
g_original <- map %*% cell_means[fit_cells,,drop=FALSE]
checks <- rbindlist(lapply(model_names,function(stem) {
  g <- if(stem=="moving_raw") g_original else sweep(g_original,1,rowMeans(g_original),"-")
  delta <- max(abs(g-GX[[stem]]))
  fit <- prcomp(g,center=TRUE,scale.=FALSE)
  data.table(check="independent_original_training_matrix_reconstruction",model=stem,
    maximum_absolute_error=delta,pass=delta<1e-10,
    PC1_loading_abs_cosine=abs(sum(fit$rotation[,1]*models[[stem]]$rotation[,1])))
}))
write_tsv(checks,"reconstruction_checks.tsv"); stopifnot(all(checks$pass))
write_tsv(cells,"bootstrap_design_cells.tsv")
message("Within-cell whole-fly bootstrap: re-estimating batch correction and both PCA bases")
set.seed(202609102)
for(b in seq_len(fly_B)) {
  bm <- t(vapply(cell_rows,function(ix) {
    sampled <- ix[sample.int(length(ix),length(ix),replace=TRUE)]
    colMeans(X[sampled,,drop=FALSE])
  },numeric(48)))
  g_raw <- map %*% bm[fit_cells,,drop=FALSE]
  for(stem in model_names) {
    g <- if(stem=="moving_raw") g_raw else sweep(g_raw,1,rowMeans(g_raw),"-")
    fit <- prcomp(g,center=TRUE,scale.=FALSE)
    z <- metrics(fit,stem,perturbed=g)
    z[,`:=`(model=stem,resampling="within_cell_whole_fly_refit_batch_and_PCA",replicate=b)]
    results[[ri<-ri+1L]] <- z
    sgn <- sign(sum(models[[stem]]$rotation[,1]*fit$rotation[,1]))
    loadings[[li<-li+1L]] <- data.table(model=stem,resampling=z$resampling,replicate=b,
      bin=1:48,loading=fit$rotation[,1]*sgn)
  }
  if(b%%50L==0L) message("  fly bootstrap ",b,"/",fly_B)
}
res <- rbindlist(results)
write_tsv(res,"pca_bootstrap_replicates.tsv")
measure_cols <- setdiff(names(res),c("model","resampling","replicate"))
numeric_res <- copy(res)
numeric_res[,(measure_cols):=lapply(.SD,as.numeric),.SDcols=measure_cols]
long <- melt(numeric_res,id.vars=c("model","resampling","replicate"),measure.vars=measure_cols,
  variable.name="metric",value.name="value")
summary <- long[,.(n_valid=sum(is.finite(value)),median=finite_median(value),
  q025=if(any(is.finite(value))) quantile(value,.025,na.rm=TRUE) else NA_real_,
  q975=if(any(is.finite(value))) quantile(value,.975,na.rm=TRUE) else NA_real_),
  by=.(model,resampling,metric)]
write_tsv(summary,"pca_bootstrap_summary.tsv")
write_tsv(res[,.(replicates=.N,PC1_axis_mixing_count=sum(PC1_best_matching_frozen_axis!=1L),
  PC1_axis_mixing_fraction=mean(PC1_best_matching_frozen_axis!=1L)),by=.(model,resampling)],"axis_mixing.tsv")
ld <- rbindlist(loadings)
write_tsv(ld[,.(zt=(bin[1]-.5)/2,median=median(loading),q025=quantile(loading,.025),
  q975=quantile(loading,.975)),by=.(model,resampling,bin)],"PC1_loading_bootstrap_intervals.tsv")

message("Conditional candidate-panel rank sensitivity: frozen QA pool and correlation matrix")
quality <- fread(file.path(src,"qa/phenotype_quality_rank.tsv"))
pool <- quality[review_status=="priority candidate" & is.finite(max_LOO_rho_change) & max_LOO_rho_change<.15]
original <- fread(file.path(src,"candidate_panel.tsv"))
ct <- fread(file.path(src,"qa/phenotype_spearman_matrix.tsv"))
cm <- as.matrix(ct[,-1]); rownames(cm)<-ct$phenotype
primary <- c("DarkSleep__D6_D1","TotalSleep__State_D2_D6","PCA_moving_centered_shape_PC1",
  "PWakeDark__State_D2_D6","PDozeDark__State_D2_D6","TwoHarmonicR2__State_D2_D6")
greedy <- function(order_names) {
  selected <- "DarkSleep__D6_D1"
  fc <- setNames(rep(0L,6),c("sleep_architecture","LD_profile","LD_transition",
    "sleep_amount_allocation","profile_PCA","activity")); fc["sleep_amount_allocation"]<-1L
  for(pn in order_names) {
    if(pn %chin% selected || length(selected)>=10L || !pn %in% colnames(cm)) next
    family <- quality[phenotype==pn,family]
    if(family %in% names(fc) && fc[family]>=2L) next
    if(any(abs(cm[pn,intersect(selected,colnames(cm))])>.85,na.rm=TRUE)) next
    selected <- c(selected,pn); if(family %in% names(fc))fc[family]<-fc[family]+1L
  }
  selected
}
stopifnot(identical(greedy(pool$phenotype),original$phenotype))
pairs <- fread(file.path(src,"qa/repeated_line_pairs.tsv"))
gp <- sort(unique(pairs$Genotype))
audit_traits <- unique(c(pool$phenotype,primary))
pp <- lapply(setNames(audit_traits,audit_traits),function(pn) {
  z <- pairs[phenotype==pn]; list(x=z$x,y=z$y,ix=split(seq_len(nrow(z)),z$Genotype))
})
set.seed(202609103)
selection <- list(); ranks <- list(); boot_rows <- list()
for(b in seq_len(selection_B)) {
  draw <- sample(gp,length(gp),replace=TRUE)
  rho <- vapply(pp,function(p) {
    ix <- unlist(p$ix[draw],use.names=FALSE)
    safe_cor(c(p$x[ix],p$y[ix]),c(p$y[ix],p$x[ix]),"spearman")
  },numeric(1))
  ordered <- pool$phenotype[order(-rho[pool$phenotype],seq_len(nrow(pool)),na.last=TRUE)]
  ordered <- ordered[is.finite(rho[ordered])]
  chosen <- greedy(ordered)
  selection[[b]] <- data.table(replicate=b,selection_position=seq_along(chosen),phenotype=chosen)
  ranks[[b]] <- data.table(replicate=b,phenotype=pool$phenotype,
    rank=match(pool$phenotype,ordered))
  boot_rows[[b]] <- data.table(replicate=b,phenotype=names(rho),repeatability_rho=rho)
}
sel <- rbindlist(selection); rk <- rbindlist(ranks); rb <- rbindlist(boot_rows)
write_tsv(sel,"conditional_panel_selections.tsv")
write_tsv(rb,"conditional_repeatability_bootstrap.tsv")
family_selection <- merge(sel,quality[,.(phenotype,family)],by="phenotype")
family_grid <- CJ(replicate=seq_len(selection_B),family=unique(quality$family))
family_selection <- merge(family_grid,family_selection[,.N,by=.(replicate,family)],
  by=c("replicate","family"),all.x=TRUE)
family_selection[is.na(N),N:=0L]
write_tsv(family_selection[,.(fraction_represented=mean(N>0),median_number_selected=median(N),
  minimum_number_selected=min(N),maximum_number_selected=max(N)),by=family],
  "conditional_family_representation.tsv")
ss <- merge(quality[phenotype %chin% audit_traits,.(phenotype,family,repeatability_rho)],
  sel[,.(selected_replicates=.N,selection_frequency=.N/selection_B),by=phenotype],by="phenotype",all.x=TRUE)
ss[is.na(selected_replicates),`:=`(selected_replicates=0L,selection_frequency=0)]
ss <- merge(ss,rk[,.(rank_median=median(rank,na.rm=TRUE),rank_q025=quantile(rank,.025,na.rm=TRUE),
  rank_q975=quantile(rank,.975,na.rm=TRUE)),by=phenotype],by="phenotype",all.x=TRUE)
ss <- merge(ss,rb[,.(rho_q025=quantile(repeatability_rho,.025,na.rm=TRUE),
  rho_q975=quantile(repeatability_rho,.975,na.rm=TRUE)),by=phenotype],by="phenotype",all.x=TRUE)
ss[,`:=`(in_frozen_QA_pool=phenotype %chin% pool$phenotype,
  in_original_discovery_panel=phenotype %chin% original$phenotype,
  frozen_primary_GWAS=phenotype %chin% primary)]
setorder(ss,-selection_frequency,-repeatability_rho)
write_tsv(ss,"conditional_selection_summary.tsv")
write_tsv(ss[frozen_primary_GWAS==TRUE],"frozen_primary_trait_selection_context.tsv")
overlap <- sel[,.(original_panel_overlap=sum(phenotype %chin% original$phenotype),panel_size=.N),by=replicate]
write_tsv(overlap,"conditional_original_panel_overlap.tsv")

fmt <- function(stem,scheme,metric_name) {
  z<-summary[model==stem & resampling==scheme & as.character(metric)==metric_name]
  sprintf("%.3f [%.3f, %.3f]",z$median,z$q025,z$q975)
}
lines <- c("# PCA and candidate-panel stability sensitivity", "",
  "This module preserves the original six primary GWAS traits, phenotype definitions, QC cohort and association results. No SNPs were used in any refit.","",
  "## Completed analyses", "",
  sprintf("- %d DGRP-line bootstrap refits for each of moving_raw and moving_centered_shape. Batch-adjusted training profiles are frozen; lines are sampled with replacement and equal draw weights. Out-of-bag correlations project held-out line profiles but their batch adjustment was fitted previously, so this is NOT independent predictive validation.",line_B),
  sprintf("- %d within-line-batch whole-fly bootstrap refits for both bases. Whole 48-bin D2-D6 profiles, not bins or days, are sampled. The complete additive Batch + Genotype multivariate model and equal-line PCA are refitted. The observed cells, fly counts, eligibility and batch design remain fixed. This is conditional sampling stability, not a new-batch/monitor or QC-uncertainty test.",fly_B),
  "- A separately constructed multivariate EMM projection reproduces both frozen training matrices to <1e-10 before resampling (reconstruction_checks.tsv). This validates the optimized recalibration against the original per-bin implementation.",
  sprintf("- %d jointly resampled genotype-cluster repeatability/ranking runs reapply the original greedy discovery-panel rule. The original QA pool (including its CI and LOO gates), phenotype values/PCA bases, correlation matrix, two-per-family cap, ten-trait cap and forced night-change anchor are held fixed. These are conditional selection frequencies, NOT full feature-selection probabilities or an independent validation of the six primary traits.",selection_B),"",
  "## PCA results", "", "Values below are bootstrap medians [2.5th, 97.5th percentile]; these summarize perturbation stability, not biological effect confidence intervals.", "")
for(stem in model_names) for(scheme in unique(res$resampling)) lines <- c(lines,
  sprintf("- %s / %s: PC1 loading cosine %s; fixed-profile PC1 score rank concordance %s; first-two-PC minimum subspace cosine %s.",stem,scheme,
    fmt(stem,scheme,"PC1_loading_cosine"),fmt(stem,scheme,"PC1_fixed_profile_score_rho"),fmt(stem,scheme,"first2_subspace_min_cosine")))
lines <- c(lines,"",
  "Fixed-profile score concordance isolates the effect of changing the PCA basis. The following fly-bootstrap score correlations instead include both changed line estimates/batch calibration and the changed PCA basis:","",
  vapply(model_names,function(stem) sprintf("- %s: fully perturbed PC1 line-score rank concordance %s.",
    stem,fmt(stem,"within_cell_whole_fly_refit_batch_and_PCA","PC1_perturbed_profile_score_rho")),character(1)))
lines <- c(lines,"", "## Interpreting the panel audit", "",
  sprintf("The frozen discovery panel has %d traits; the frozen QA pool has %d. Median original-panel overlap across rank perturbations is %.1f/%d (range %d-%d).",nrow(original),nrow(pool),median(overlap$original_panel_overlap),nrow(original),min(overlap$original_panel_overlap),max(overlap$original_panel_overlap)),
  "The manually consolidated six-primary-trait set is not identical to the discovery panel: TotalSleep replaces a highly redundant raw-moving PC1 while centered-shape PC1 represents timing/profile structure. Low greedy inclusion is therefore not a rationale to remove a preselected primary trait; high inclusion is not biological validation. The night-change anchor has selection frequency one by construction.",
  "Read conditional_selection_summary.tsv alongside the frozen primary-trait context table. Rank exchanges between correlated or same-family alternatives can change membership without contradicting evidence for the broader phenotype dimension.","",
  "## Limits and next evidence", "",
  "Centered shape subtracts each profile mean but does not normalize amplitude. Good loading stability does not establish an endogenous circadian oscillator, causal isolation response, or a clock knockout-like mechanism. LD rhythm descriptors can include masking and transitions.",
  "Fly bootstrap treats flies as sampling units conditional on observed line-batch cells; it cannot resolve line-monitor confounding or estimate independent monitor/batch perturbations. Neither bootstrap includes death-rule rescue or new flies/lines; that is a separate QC/replication question.",
  "Primary confirmatory inference must not use these data-driven discovery selection frequencies as preregistration or multiplicity correction. No new GWAS tests or p-value thresholds were introduced.","",
  "## Reproduction", "", "Run from the project root:", "", "```sh", "Rscript scripts/readiness_pca.R", "```", "",
  "Input checksums, deterministic seeds, complete replicate tables and sessionInfo.txt are included. Original files are read-only inputs.")
writeLines(lines,file.path(out,"README.md"))
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
validation <- data.table(check=c("input_checksums_unchanged","expected_PCA_replicate_count",
  "unit_loading_cosines_in_range","original_greedy_panel_exactly_reproduced",
  "all_selection_replicates_complete","independent_calibration_matches_frozen"),
  pass=c(identical(unname(tools::md5sum(input_files)),manifest$md5),
    nrow(res)==2L*(line_B+fly_B),
    all(res$PC1_loading_cosine<=1+1e-10 & res$PC1_loading_cosine>=0),
    identical(greedy(pool$phenotype),original$phenotype),
    uniqueN(sel$replicate)==selection_B && all(overlap$panel_size==10L),all(checks$pass)))
write_tsv(validation,"validation_checks.tsv")
stopifnot(all(validation$pass))
message("Completed PCA stability and conditional panel audit: ",out)
