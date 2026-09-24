#!/usr/bin/env Rscript
# Run from the release root. Preserves effects and applies the published mask.
suppressPackageStartupMessages(library(data.table))
setDTthreads(2)
root <- 'reproducibility/profile_gwas_call90_maf10'
out <- 'rerun/eligible_gwas';dir.create(out,recursive=TRUE,showWarnings=FALSE)
plan <- fread(file.path(root,'inputs/run_plan.tsv'))
expected <- fread(file.path(root,'results/gwas_summary.tsv'))
get_counts <- function(path,n) {
  z <- fread(path)
  setnames(z,c('chr','SNP','A1','A2','hom1','het','hom2','hap1','hap2','missing'))
  stopifnot(all(z$hom1+z$het+z$hom2+z$missing==n),all(z$hap1+z$hap2==0))
  z[,called:=hom1+het+hom2]
  z[,frequency:=(hom1+het/2)/called]
  z[,`:=`(MAF=pmin(frequency,1-frequency),call_rate=called/n)]
  z[is.finite(MAF) & MAF>=.1-1e-12 & call_rate>=.9-1e-12,.(SNP,MAF,call_rate)]
}
common <- get_counts(file.path(root,'shared/cohort_counts.frqx'),176)
checks <- list()
for(i in seq_len(nrow(plan))) {
  z <- plan[i];stem <- paste(z$trait,z$variant,sep='__')
  counts <- if(z$variant=='exclude_influential')
    get_counts(file.path(root,'runs',paste0(stem,'.counts.frqx')),z$n_lines) else common
  extension <- if(z$mode=='loco') '.loco.mlma' else '.mlma'
  a <- fread(file.path(root,'runs',paste0(stem,extension)))
  a <- merge(a,counts,by='SNP',sort=FALSE)
  sm <- expected[trait==z$trait & variant==z$variant]
  stopifnot(nrow(a)==sm$eligible_variants,all(is.finite(a$p)),all(a$p>0 & a$p<=1),
            sum(a$p<1e-5)==sm$nominal_hits,min(a$p)==sm$minimum_p)
  inflation <- median(qchisq(a$p,1,lower.tail=FALSE))/qchisq(.5,1)
  stopifnot(abs(inflation-sm$lambda_gc)<1e-10)
  fwrite(a,file.path(out,paste0(stem,'.tsv.gz')),sep='\t',compress='gzip')
  checks[[i]] <- data.table(trait=z$trait,variant=z$variant,eligible_variants=nrow(a),
                             nominal_hits=sum(a$p<1e-5),minimum_p=min(a$p),lambda_gc=inflation)
  message('Verified ',stem)
}
fwrite(rbindlist(checks),file.path(out,'verification.tsv'),sep='\t')
