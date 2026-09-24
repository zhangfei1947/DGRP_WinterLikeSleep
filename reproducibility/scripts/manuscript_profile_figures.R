#!/usr/bin/env Rscript
# Presentation-only rebuild from frozen derived data. Run from project root.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(ggrepel)})
setDTthreads(2); set.seed(20260911)
out <- 'manuscript_figures_tables'
for(d in c('figures/main','figures/supplementary','source_data','qa')) dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
inputs <- character()
read_tab <- function(p) {inputs <<- union(inputs,p); fread(p)}
read_rds <- function(p) {inputs <<- union(inputs,p); as.data.table(readRDS(p))}
write_data <- function(x,n) fwrite(x,file.path(out,'source_data',paste0(n,'.tsv')),sep='\t',na='NA')
blue <- '#376B8C'; orange <- '#B66B36'; ink <- '#252525'
theme_set(theme_classic(base_size=10,base_family='Helvetica')+theme(
 plot.title=element_text(size=11,face='bold',hjust=0),plot.subtitle=element_text(size=9),
 plot.tag=element_text(size=12,face='bold'),plot.caption=element_text(size=8,hjust=0),
 axis.text=element_text(colour=ink,size=9),axis.title=element_text(size=10),
 strip.background=element_blank(),strip.text=element_text(face='bold',size=10),
 legend.position='top',legend.title=element_blank(),legend.text=element_text(size=9),
 panel.spacing.x=grid::unit(.28,'inches'),plot.margin=margin(7,9,7,7)))
save_plot <- function(p,n,w=7.2,h=6,supp=FALSE) {
  stem <- file.path(out,'figures',if(supp)'supplementary' else 'main',n)
  ggsave(paste0(stem,'.pdf'),p,width=w,height=h,device=pdf,useDingbats=FALSE)
  ggsave(paste0(stem,'.png'),p,width=w,height=h,dpi=300,bg='white')
}
ldbar <- function(ymin=-.09,ymax=-.045) list(
 annotate('rect',xmin=0,xmax=8,ymin=ymin,ymax=ymax,fill='white',colour='black',linewidth=.3),
 annotate('rect',xmin=8,xmax=24,ymin=ymin,ymax=ymax,fill='black',colour='black',linewidth=.3))
traits <- read_tab('profile_gwas_call90_maf10/inputs/primary_traits.tsv')
labels <- c(NightSleepChange='Night sleep change',MovingShapePC1='Moving shape PC1',
 NightPWake='Night P(Wake)',NightPDoze='Night P(Doze)',LDHarmonicR2='LD harmonic R2',TotalSleep='Total sleep')
traits[,short:=unname(labels[trait])]
quality <- read_tab('profile_phenotype_discovery/qa/phenotype_quality_rank.tsv')
line <- read_tab('profile_phenotype_discovery/data/line_phenotypes_batchFE.tsv')
lb <- read_tab('profile_phenotype_discovery/data/line_batch_phenotypes.tsv')
inventory <- read_tab('manuscript_readiness/provenance/batch_inventory.tsv')
daily <- read_rds('profile_phenotype_discovery/data/individual_daily_features.rds')
pop <- read_tab('profile_phenotype_discovery/data/population_profiles_D2_D6.tsv')
pairs <- read_tab('profile_phenotype_discovery/qa/repeated_line_pairs.tsv')
loadings <- read_tab('profile_phenotype_discovery/pca/loadings.tsv')
stopifnot(nrow(traits)==6,all(traits$phenotype %in% line$phenotype),sum(inventory$retained)==uniqueN(daily$id))

# Figure 1: the design is not an isolation-only contrast. Daily summaries use
# identical flies on all six days within each phase; no D1 imputation.
timeline <- ggplot()+
 annotate('segment',x=.4,xend=6.6,y=2,yend=2,linewidth=.5,arrow=grid::arrow(length=grid::unit(.09,'inches')))+
 annotate('point',x=c(1.1,3.4,5.8),y=2,size=2)+
 annotate('text',x=1.1,y=2.65,label='Group housing\n2-5 days after eclosion',size=3.3,lineheight=1.1)+
 annotate('text',x=1.1,y=1.3,label='12L:12D\n21.5 C',size=3.3)+
 annotate('text',x=3.4,y=2.65,label='D0: male sorting\nand individual DAM loading',size=3.3,lineheight=1.1)+
 annotate('text',x=3.4,y=1.3,label='Group to individual\n12L:12D to 8L:16D',size=3.3)+
 annotate('text',x=5.8,y=2.65,label='D1-D6 recording\nIndividual housing',size=3.3,lineheight=1.1)+
 annotate('text',x=5.8,y=1.3,label='8L:16D\n21.5 C',size=3.3)+
 coord_cartesian(xlim=c(0,7),ylim=c(.7,3.1),clip='off')+theme_void(base_family='Helvetica')+
 labs(title='Housing and light-cycle transition')
phases <- c(DarkSleep='Night sleep (ZT8-24)',LightSleep='Day sleep (ZT0-8)')
daily_cells <- rbindlist(lapply(names(phases),function(m){
 d<-daily[,.(id,Batch,Genotype,day,value=get(m))]
 ids<-d[,.(complete=.N==6 && all(is.finite(value))),by=id][complete==TRUE,id]
 d<-d[id%chin%ids]
 cells<-d[,.(value=median(value),n_flies=.N),by=.(Batch,Genotype,day)][n_flies>=8]
 cells[,phase:=phases[[m]]]; cells
}))
daily_lines<-daily_cells[Genotype!='CS',.(value=mean(value),n_batches=.N),by=.(Genotype,day,phase)]
daily_summary<-daily_lines[,.(median=median(value),q25=quantile(value,.25),q75=quantile(value,.75),n_lines=.N),by=.(day,phase)]
daily_summary[,phase:=factor(phase,levels=unname(phases))]
dp<-ggplot(daily_summary,aes(day,median))+geom_ribbon(aes(ymin=q25,ymax=q75),fill=blue,alpha=.15)+
 geom_line(colour=blue,linewidth=.7)+geom_point(colour=blue,size=1.8)+
 facet_wrap(~phase,ncol=2,scales='free_y',labeller=as_labeller(c('Night sleep (ZT8-24)'='Night sleep (ZT8-24), n = 176','Day sleep (ZT0-8)'='Day sleep (ZT0-8), n = 169')))+scale_x_continuous(breaks=1:6,labels=paste0('D',1:6))+
 labs(x=NULL,y='Sleep (min)',title='Daily sleep across DGRP lines',subtitle='Median and interquartile range of line summaries')
counts<-data.table(stage=c('Planned fly records','Retained flies','Eligible DGRP flies','Eligible DGRP lines'),
 n=c(sum(inventory$planned),sum(inventory$retained),sum(lb[phenotype=='DarkSleep__D6_D1' & Genotype!='CS' & n_flies>=8,n_flies]),
 line[phenotype=='DarkSleep__D6_D1',.N]))
counts[,stage:=factor(stage,levels=rev(stage))]
cohort<-ggplot(counts,aes(n,stage))+geom_point(size=2.3,colour=blue)+geom_text(aes(label=format(n,big.mark=',')),hjust=-.25,size=3.3)+
 scale_x_continuous(limits=c(0,max(counts$n)*1.16),breaks=c(0,3000,6000,9000),labels=scales::comma)+
 labs(x='Count',y=NULL,title='Cohort',subtitle=sprintf('%d batch labels; eligibility shown for night sleep change',nrow(inventory)))
save_plot(timeline/cohort/dp+plot_layout(heights=c(1,1,1.3))+plot_annotation(tag_levels='A'),
 'Figure1_design_and_cohort',h=7.2)
write_data(counts,'Figure1_cohort');write_data(daily_cells,'Figure1_daily_line_batch');write_data(daily_lines,'Figure1_daily_line');write_data(daily_summary,'Figure1_daily_summary')

# Figure 2: line-level distribution and independent-batch descriptive agreement.
deltas<-line[phenotype%chin%c('DarkSleep__D6_D1','LightSleep__D6_D1')]
deltas[,phase:=ifelse(phenotype=='DarkSleep__D6_D1','Night sleep','Day sleep')]
deltas[,phase:=factor(phase,levels=c('Night sleep','Day sleep'))]
ann<-deltas[,.(label=paste0('n = ',.N,' lines')),by=phase]
hist<-ggplot(deltas,aes(value))+geom_histogram(binwidth=40,boundary=0,fill=blue,colour='white',linewidth=.25)+
 geom_vline(xintercept=0,linetype='dashed',colour='grey40',linewidth=.4)+
 geom_text(data=ann,aes(x=Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=1.05,vjust=1.4,size=3.1)+
 facet_wrap(~phase,ncol=2)+labs(x='Batch-adjusted D6 - D1 sleep (min)',y='DGRP lines',title='Variation in sleep change')
reps<-pairs[phenotype%chin%deltas$phenotype]
repplots<-lapply(c('DarkSleep__D6_D1','LightSleep__D6_D1'),function(pn){
 x<-reps[phenotype==pn];q<-quality[phenotype==pn];lims<-range(c(x$x,x$y))
 ggplot(x,aes(x,y))+geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey55',linewidth=.4)+
 geom_point(colour=blue,size=1.7,alpha=.65)+coord_equal(xlim=lims,ylim=lims)+
 labs(title=if(pn=='DarkSleep__D6_D1')'Night sleep' else 'Day sleep',
 subtitle=sprintf('%d lines; %d batch pairs; rho = %.2f',q$repeated_lines,q$batch_pairs,q$repeatability_rho),
 x='First batch in pair (min)',y='Second batch in pair (min)')
})
save_plot(hist/(repplots[[1]]|repplots[[2]])+plot_layout(heights=c(.85,1))+plot_annotation(tag_levels='A'),
 'Figure2_sleep_change',h=6.8)
write_data(deltas,'Figure2_line_sleep_change');write_data(reps,'Figure2_repeated_batch_pairs')

# Figure 3: observed probabilities, not PCA reconstructions. Preserve physical
# range 0-1. Equal eligible batches per line; sorting only uses frozen shape PC1.
scores<-line[phenotype=='PCA_moving_centered_shape_PC1',.(Genotype,PC1=value)]
setorder(scores,PC1,Genotype);scores[,rank:=.I]
pp<-pop[n_flies>=8 & Genotype%chin%scores$Genotype]
atlas<-pp[,.(value=mean(value),n_batches=.N,total_flies=sum(n_flies)),by=.(Genotype,slot,signal)]
atlas<-merge(atlas,scores,by='Genotype');atlas[,zt:=(slot+.5)/2]
stopifnot(nrow(atlas)==176*48*2,all(atlas$value>=0 & atlas$value<=1))
heat<-lapply(c('moving','sleep'),function(sig){
 ggplot(atlas[signal==sig],aes(zt,rank,fill=value))+geom_tile(width=.5,height=1)+
 scale_fill_gradient(low='white',high=blue,limits=c(0,1),breaks=c(0,.5,1),name='Probability')+
 annotate('rect',xmin=0,xmax=8,ymin=-8,ymax=-3,fill='white',colour='black',linewidth=.3)+
 annotate('rect',xmin=8,xmax=24,ymin=-8,ymax=-3,fill='black',colour='black',linewidth=.3)+
 scale_x_continuous(breaks=c(0,8,16,24),expand=c(0,0))+
 scale_y_continuous(breaks=c(1,88,176),limits=c(-8,176.5),expand=c(0,0))+
 labs(title=if(sig=='moving')'Moving profiles' else 'Sleep profiles',x='ZT (h)',y='DGRP line rank (shape PC1)')+
 theme(legend.position='bottom',legend.title=element_text(size=9),axis.line=element_blank(),axis.ticks.y=element_blank())
})
# Display selection rule: ranks nearest 10%, 50%, 90% of frozen PC1 order.
examples<-scores[round(c(.1,.5,.9)*(nrow(scores)-1))+1]
examples[,selection:=c('PC1 10th percentile','PC1 median','PC1 90th percentile')]
ep<-pop[signal=='moving' & n_flies>=8 & Genotype%chin%c('CS',examples$Genotype)]
ep[,zt:=(slot+.5)/2]
elabs<-ep[,.(nb=uniqueN(Batch),n=sum(n_flies[slot==0])),by=Genotype]
elabs[,label:=sprintf('%s\nn = %d flies, %d %s',Genotype,n,nb,ifelse(nb==1,'batch','batches'))]
ep<-merge(ep,elabs,by='Genotype')
ep[,label:=factor(label,levels=elabs[match(c('CS',examples$Genotype),Genotype),label])]
explot<-ggplot(ep,aes(zt,value,group=Batch))+geom_line(colour=blue,alpha=.55,linewidth=.5)+ldbar()+
 facet_wrap(~label,ncol=2)+scale_x_continuous(limits=c(0,24),breaks=c(0,8,16,24),expand=c(.01,0))+
 scale_y_continuous(limits=c(-.095,1),breaks=c(0,.5,1),expand=c(0,0))+
 labs(title='Observed D2-D6 moving profiles',subtitle='One curve per eligible batch; CS and PC1-selected DGRP examples',x='ZT (h)',y='Moving probability')+
 theme(strip.text=element_text(size=9))
save_plot((heat[[1]]|heat[[2]])/explot+plot_layout(heights=c(1.1,1.3))+plot_annotation(tag_levels='A'),
 'Figure3_LD_profile_diversity',h=8.4)
write_data(atlas,'Figure3_profile_atlas');write_data(examples,'Figure3_example_selection');write_data(ep,'Figure3_example_batch_profiles')

# Figure 4: six primary phenotypes and conditional PCA interpretation.
primary<-merge(traits,line,by='phenotype');stopifnot(primary[,.N,by=trait][,all(N==176)])
mat<-dcast(primary,Genotype~trait,value.var='value');M<-as.matrix(mat[,traits$trait,with=FALSE]);rownames(M)<-mat$Genotype
C<-cor(M,method='spearman');corr<-as.data.table(as.table(C));setnames(corr,c('a','b','rho'))
corr[,a:=factor(a,levels=traits$trait,labels=traits$short)];corr[,b:=factor(b,levels=rev(traits$trait),labels=rev(traits$short))]
cp<-ggplot(corr,aes(a,b,fill=rho))+geom_tile(colour='white',linewidth=.5)+geom_text(aes(label=sprintf('%.2f',rho)),size=2.8)+
 scale_fill_gradient2(low=orange,mid='white',high=blue,limits=c(-1,1),name='Spearman rho')+
 coord_equal()+labs(x=NULL,y=NULL,title='Primary phenotype correlations')+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.text.y=element_text(size=8),legend.position='bottom')
rq<-merge(traits,quality,by='phenotype');rq[,short:=factor(short,levels=rev(traits$short))]
rp<-ggplot(rq,aes(repeatability_rho,short))+geom_segment(aes(x=rho_ci_low,xend=rho_ci_high,yend=short),linewidth=.55,colour='grey55')+
 geom_point(colour=blue,size=2.3)+scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1))+
 labs(x='Repeatability (Spearman rho)',y=NULL,title='Cross-batch agreement',subtitle='95% cluster-bootstrap intervals')
ld<-loadings[model=='moving_centered_shape' & PC=='PC1']
yr<-range(ld$loading);dy<-diff(yr)
lp<-ggplot(ld,aes(zt,loading))+geom_hline(yintercept=0,colour='grey70',linewidth=.3)+geom_line(colour=blue,linewidth=.7)+
 ldbar(yr[1]-.15*dy,yr[1]-.08*dy)+scale_x_continuous(limits=c(0,24),breaks=c(0,8,16,24))+
 labs(title='Moving shape PC1',subtitle=sprintf('%.1f%% of profile variance',100*ld$variance_fraction[1]),x='ZT (h)',y='PC1 loading')
sc<-mat[,.(Genotype,MovingShapePC1,NightSleepChange)]
rho<-cor(sc$MovingShapePC1,sc$NightSleepChange,method='spearman')
sp<-ggplot(sc,aes(MovingShapePC1,NightSleepChange))+geom_point(colour=blue,alpha=.65,size=1.4)+
 geom_hline(yintercept=0,colour='grey65',linetype='dashed',linewidth=.3)+
 labs(title='Shape versus sleep change',subtitle=sprintf('n = %d lines; rho = %.2f',nrow(sc),rho),x='Moving shape PC1',y='Night sleep D6 - D1 (min)')
save_plot((cp|rp)/(lp|sp)+plot_layout(heights=c(1.3,1))+plot_annotation(tag_levels='A'),'Figure4_primary_phenotypes',h=7.8)
write_data(corr,'Figure4_correlations');write_data(rq,'Figure4_repeatability');write_data(ld,'Figure4_PC1_loadings');write_data(sc,'Figure4_shape_vs_change')

# Supplementary S1: controls are internal descriptive controls, not an
# independently replicated isolation experiment.
cs<-pop[Genotype=='CS' & n_flies>=8];cs[,zt:=(slot+.5)/2]
cs[,signal:=factor(signal,levels=c('moving','sleep'),labels=c('Moving','Sleep'))]
csp<-ggplot(cs,aes(zt,value,group=Batch))+geom_line(colour=blue,alpha=.35,linewidth=.5)+ldbar()+
 facet_wrap(~signal,ncol=1)+scale_x_continuous(limits=c(0,24),breaks=c(0,8,16,24))+
 scale_y_continuous(limits=c(-.095,1),breaks=c(0,.5,1))+
 labs(title='CS population profiles across batches',subtitle=sprintf('D2-D6; %d eligible batches; one curve per batch',uniqueN(cs$Batch)),x='ZT (h)',y='Probability')
save_plot(csp,'FigureS1_CS_profiles',h=6.2,supp=TRUE);write_data(cs,'FigureS1_CS_profiles')

a<-read_tab('manuscript_readiness/baseline/time_window_common_fly_cohort_agreement.tsv')
a[,trait:=factor(metric,levels=c('DarkSleep','LightSleep','TotalSleep'),labels=c('Night sleep','Day sleep','Total sleep'))]
a[,window:=factor(contrast,levels=c('Late2_D2D3','D6_D2','Late2_Early2'),labels=c('D5-D6 mean minus D2-D3 mean','D6 minus D2','D5-D6 mean minus D1-D2 mean'))]
tp<-ggplot(a,aes(spearman,window))+geom_point(colour=blue,size=2.2)+
 geom_text(aes(label=sprintf('%.3f (n=%d)',spearman,n_common_lines)),hjust=-.15,size=3)+facet_wrap(~trait,ncol=1)+
 scale_x_continuous(limits=c(0,1.2),breaks=c(0,.25,.5,.75,1))+
 labs(title='Sleep-change rankings across time windows',subtitle='Each alternative versus D6-D1 in matched fly cohorts',x='Spearman correlation of line estimates',y=NULL)
save_plot(tp,'FigureS2_time_windows',h=6.2,supp=TRUE);write_data(a,'FigureS2_time_windows')

qc<-read_tab('manuscript_readiness/qc/concordance.tsv')[estimate=='batchFE']
qc<-merge(qc,traits[,.(phenotype,short)],by='phenotype')
scheme_labels<-c(current_simulated='Reimplemented QC',recovery_aware='Recovery-aware',prop005_24h='Activity fraction 0.005 / 24 h',prop001_24h='Activity fraction 0.001 / 24 h',zero_24h='Zero movement / 24 h',prop01_48h='Activity fraction 0.01 / 48 h')
stopifnot(all(qc$scheme%in%names(scheme_labels)))
qc[,rule:=factor(scheme,levels=rev(names(scheme_labels)),labels=rev(scheme_labels))]
qc[,short:=factor(short,levels=traits$short)]
qp<-ggplot(qc,aes(rho,rule))+geom_point(colour=blue,size=2)+
 facet_wrap(~short,ncol=2)+scale_x_continuous(limits=c(.95,1),breaks=c(.95,.975,1))+
 labs(title='Primary phenotype rankings after QC changes',subtitle='Compared with recorded QC; common lines; frozen PC basis',x='Spearman correlation (focused scale)',y=NULL)
save_plot(qp,'FigureS3_QC_sensitivity',h=7.2,supp=TRUE);write_data(qc,'FigureS3_QC_concordance')

boot<-read_tab('manuscript_readiness/pca/PC1_loading_bootstrap_intervals.tsv')
ref<-loadings[PC=='PC1' & model%chin%unique(boot$model),.(model,zt,reference=loading)]
boot<-merge(boot,ref,by=c('model','zt'))
boot[,resampling_label:=ifelse(grepl('DGRP_line',resampling),'Line resampling','Whole-fly resampling')]
boot[,model_label:=ifelse(model=='moving_raw','Raw moving PC1','Centered moving PC1')]
strips<-boot[,.(ymin=min(q025)-.09*diff(range(c(q025,q975))),ymax=min(q025)-.04*diff(range(c(q025,q975)))),by=model_label]
strips<-merge(strips,unique(boot[,.(model_label,resampling_label)]),by='model_label')
bp<-ggplot(boot,aes(zt,median))+geom_ribbon(aes(ymin=q025,ymax=q975),fill=blue,alpha=.17)+
 geom_line(colour=blue,linewidth=.6)+geom_line(aes(y=reference),linetype='dashed',colour='grey35',linewidth=.4)+
 geom_hline(data=unique(boot[model=='moving_centered_shape',.(model_label,resampling_label)])[,zero:=0],aes(yintercept=zero),colour='grey75',linewidth=.3)+
 geom_rect(data=strips,aes(xmin=0,xmax=8,ymin=ymin,ymax=ymax),inherit.aes=FALSE,fill='white',colour='black',linewidth=.3)+
 geom_rect(data=strips,aes(xmin=8,xmax=24,ymin=ymin,ymax=ymax),inherit.aes=FALSE,fill='black')+
 facet_grid(model_label~resampling_label,scales='free_y')+scale_x_continuous(limits=c(0,24),breaks=c(0,8,16,24))+
 labs(title='PC1 loading stability',subtitle='Solid: bootstrap median; dashed: frozen loading',x='ZT (h)',y='PC1 loading',caption='Bands: pointwise 2.5th-97.5th percentiles; not simultaneous confidence bands.')
save_plot(bp,'FigureS4_PCA_stability',h=6,supp=TRUE);write_data(boot,'FigureS4_PCA_stability')

sel<-read_tab('manuscript_readiness/pca/conditional_selection_summary.tsv')[in_original_discovery_panel==TRUE]
shorts<-c(DarkSleep__D6_D1='Night sleep change (forced)',PCA_moving_centered_shape_PC1='Moving shape PC1',PCA_sleep_centered_shape_PC1='Sleep shape PC1',
 TwoHarmonicR2__State_D2_D6='LD harmonic R2',PDozeDark__State_D2_D6='Night P(Doze)',ActivityLight__State_D2_D6='Day activity rate',
 PCA_moving_raw_PC1='Raw moving PC1',LightSleep__State_D2_D6='Day sleep',Harmonic24Amplitude__State_D2_D6='24-h harmonic amplitude',
 PWakeDark__State_D2_D6='Night P(Wake)',MovingLight__D6_D1='Day moving D6 - D1')
sel[,label:=gsub('__State_D2_D6',' (D2-D6)',phenotype,fixed=TRUE)]
sel[phenotype%in%names(shorts),label:=unname(shorts[phenotype])]
sel[,label:=factor(label,levels=label[order(selection_frequency)])]
ss<-ggplot(sel,aes(selection_frequency,label))+geom_point(colour=blue,size=2.4)+
 scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.25),labels=scales::percent)+
 labs(title='Conditional discovery-panel membership',subtitle='500 genotype-cluster rank resamples; original 10-trait panel',x='Inclusion frequency',y=NULL,
 caption='The night sleep anchor was forced. This is not selection probability for the final six primary traits.')
save_plot(ss,'FigureS6_discovery_panel_stability',h=5.5,supp=TRUE);write_data(sel,'FigureS6_discovery_panel_stability')

fwrite(data.table(source=inputs,md5=unname(tools::md5sum(inputs))),file.path(out,'qa','profile_input_manifest.tsv'),sep='\t')
checks<-data.table(check=c('retained_flies_match_inventory','six_primary_traits_each_176_lines','atlas_176_lines_48bins_2signals','atlas_probability_bounds','PC1_shape_sleep_change_rho'),
 value=c(as.character(uniqueN(daily$id)),'176',as.character(nrow(atlas)),'0 to 1',sprintf('%.8f',rho)))
fwrite(checks,file.path(out,'qa','profile_checks.tsv'),sep='\t')
capture.output(sessionInfo(),file=file.path(out,'qa','profile_sessionInfo.txt'))
message('Built four main and five supplementary figures in ',out)
