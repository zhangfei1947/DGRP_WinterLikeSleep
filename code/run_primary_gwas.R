#!/usr/bin/env Rscript
# Run from release root. Set GCTA_BIN to the GCTA 1.95.3 executable.
suppressPackageStartupMessages(library(data.table))
setDTthreads(2)
exe <- Sys.getenv('GCTA_BIN', 'gcta64')
args <- commandArgs(TRUE)
selected_variant <- if(length(args)) args[1] else 'primary'
allowed <- c('primary','raw_batch','RINT','exclude_influential','legacy_residual','chromosome_LOCO')
stopifnot(selected_variant %in% allowed)
root <- normalizePath('reproducibility',mustWork=TRUE)
base <- file.path(root,'publication_analysis/gwas_shared/common_min8')
gwas <- file.path(root,'profile_gwas_call90_maf10')
plan <- fread(file.path(gwas,'inputs/run_plan.tsv'))
plan <- plan[variant==selected_variant]
if(length(args)>1) plan <- plan[trait==args[2]]
stopifnot(nrow(plan)>0)
out <- file.path(Sys.getenv('GWAS_OUTPUT_DIR',file.path('rerun','gwas')),selected_variant)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(i in seq_len(nrow(plan))) {
  z <- plan[i]
  prefix <- file.path(out,paste(z$trait,z$variant,sep='__'))
  result <- paste0(prefix,if(z$mode=='loco') '.loco.mlma' else '.mlma')
  if(file.exists(result)) stop('Output already exists; use a new rerun directory: ',result)
  stem <- file.path(root,z$stem)
  bfile <- file.path(base,'plink.gwas')
  if(z$mode=='loco') {
    bfile <- file.path(out,'whole_chromosome')
    for(ext in c('bed','fam')) if(!file.exists(paste0(bfile,'.',ext)))
      stopifnot(file.copy(file.path(base,paste0('plink.gwas.',ext)),paste0(bfile,'.',ext)))
    if(!file.exists(paste0(bfile,'.bim')))
      stopifnot(file.copy(file.path(gwas,'shared/whole_chromosome.bim'),paste0(bfile,'.bim')))
  }
  a <- c(if(z$mode=='loco') '--mlma-loco' else '--mlma', '--bfile',bfile,
         '--autosome-num','32','--pheno',paste0(stem,'.pheno'),'--keep',paste0(stem,'.keep'),
         '--thread-num',Sys.getenv('GWAS_THREADS','4'),'--out',prefix)
  if(z$mode!='loco') a <- c(a,'--grm-gz',file.path(base,'grm.scaled'))
  if(z$joint_covariates) a <- c(a,'--qcovar',paste0(stem,'.qcovar'),'--mlma-no-preadj-covar')
  writeLines(c(exe,a),paste0(prefix,'.command.txt'))
  rc <- system2(exe,shQuote(a),stdout=paste0(prefix,'.console.log'),stderr=paste0(prefix,'.console.log'))
  stopifnot(rc==0,file.exists(result))
}
message('Completed ',nrow(plan),' scans. Apply the cohort-specific eligibility mask before reporting results; see docs/GWAS.md.')
