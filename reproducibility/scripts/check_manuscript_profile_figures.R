#!/usr/bin/env Rscript
# Independent arithmetic checks of the exported figures against frozen sources.
suppressPackageStartupMessages(library(data.table))
setDTthreads(2)
root<-'manuscript_figures_tables';s<-file.path(root,'source_data')
readout<-function(n) fread(file.path(s,paste0(n,'.tsv')))
checks<-list()
record<-function(label,ok,detail){stopifnot(isTRUE(ok));checks[[label]]<<-data.table(check=label,result='PASS',detail=detail)}
q<-fread('profile_phenotype_discovery/qa/phenotype_quality_rank.tsv')
p<-readout('Figure2_repeated_batch_pairs')
for(pn in unique(p$phenotype)){
 z<-p[phenotype==pn];rho<-cor(c(z$x,z$y),c(z$y,z$x),method='spearman')
 record(paste('pair rho',pn),abs(rho-q[phenotype==pn,repeatability_rho])<1e-12,sprintf('%.12f',rho))
 record(paste('pair keys',pn),uniqueN(z[,.(Genotype,Batch1,Batch2)])==nrow(z),as.character(nrow(z)))
}
traits<-fread('profile_gwas_call90_maf10/inputs/primary_traits.tsv')
line<-fread('profile_phenotype_discovery/data/line_phenotypes_batchFE.tsv')
L<-merge(traits,line,by='phenotype');W<-dcast(L,Genotype~trait,value.var='value')
expected<-fread('profile_gwas_call90_maf10/results/primary_phenotype_correlations.tsv')
delta<-mapply(function(a,b,r) abs(cor(W[[a]],W[[b]],method='spearman')-r),expected$trait_a,expected$trait_b,expected$spearman_rho)
record('six-trait correlations vs original GWAS inputs',max(delta)<1e-12,sprintf('max absolute difference %.3g',max(delta)))
pc<-readout('Figure3_example_batch_profiles')
bins<-as.data.table(readRDS('profile_phenotype_discovery/data/individual_profiles_30min.rds'))
bins<-bins[Genotype%in%pc$Genotype & day%in%2:6]
perfly<-bins[,.(value=if(.N==5 && all(is.finite(moving)))mean(moving) else NA_real_),by=.(id,Batch,Genotype,slot)]
calc<-perfly[,.(recomputed=mean(value,na.rm=TRUE),n_recomputed=sum(is.finite(value))),by=.(Batch,Genotype,slot)][n_recomputed>=8]
j<-merge(pc,calc,by=c('Batch','Genotype','slot'))
record('example profiles from individual 30-minute data',nrow(j)==nrow(pc) && max(abs(j$value-j$recomputed))<1e-12,sprintf('%d bin-by-batch rows',nrow(j)))
record('example profile sample counts',all(j$n_flies==j$n_recomputed),'All bins and selected line-batches')
atlas<-readout('Figure3_profile_atlas')
record('atlas complete probability matrix',uniqueN(atlas$Genotype)==176 && nrow(atlas)==176*48*2 && all(atlas$value>=0 & atlas$value<=1),'176 lines x 48 bins x 2 signals')
daily<-as.data.table(readRDS('profile_phenotype_discovery/data/individual_daily_features.rds'))
dc<-readout('Figure1_daily_line_batch')
# Independently check one line-batch-day median per phase.
for(m in c('DarkSleep','LightSleep')){
 phase_label<-if(m=='DarkSleep')'Night sleep (ZT8-24)' else 'Day sleep (ZT0-8)'
 one<-dc[phase==phase_label & day%in%c(1,6)][1]
 eligible<-daily[Batch==one$Batch & Genotype==one$Genotype,.(ok=.N==6 && all(is.finite(get(m)))),by=id][ok==TRUE,id]
 val<-median(daily[id%in%eligible & day==one$day,get(m)])
 record(paste('daily median sampled',m),abs(val-one$value)<1e-12,sprintf('%s %s D%d',one$Batch,one$Genotype,one$day))
}
fwrite(rbindlist(checks),file.path(root,'qa','independent_profile_checks.tsv'),sep='\t')
print(rbindlist(checks))
