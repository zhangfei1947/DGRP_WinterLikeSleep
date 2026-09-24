#!/usr/bin/env Rscript
# Render frozen v2 primary scans and fixed-lead QC sensitivity; never refit a GWAS.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
setDTthreads(2)
root <- "profile_gwas_call90_maf10"
out <- "manuscript_figures_tables"
for (d in c("figures/main","figures/supplementary","source_data"))
  dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
blue <- "#376B8C"
base_theme <- theme_classic(base_size=10,base_family="Helvetica") +
  theme(text=element_text(colour="#262626"),axis.text=element_text(colour="#333333",size=9),
        axis.title=element_text(size=10),plot.title=element_text(face="bold",size=10.5),
        plot.subtitle=element_text(size=9),plot.margin=margin(4,7,4,4),
        legend.position="none",strip.background=element_blank(),
        strip.text=element_text(face="bold",hjust=0,size=10))
theme_set(base_theme)
save_figure <- function(p,name,w,h,folder) {
  ggsave(file.path(out,"figures",folder,paste0(name,".pdf")),p,width=w,height=h,
         device=grDevices::pdf,useDingbats=FALSE,bg="white")
  ggsave(file.path(out,"figures",folder,paste0(name,".png")),p,width=w,height=h,dpi=300,bg="white")
}
write_source <- function(x,name) fwrite(x,file.path(out,"source_data",name),sep="\t",na="NA")
traits <- fread(file.path(root,"inputs/primary_traits.tsv"))
trait_order <- traits$trait
labels <- setNames(c("Night sleep change (D6 - D1)","Moving shape PC1",
                     "Night P(Wake)","Night P(Doze)","LD harmonic R2","Total sleep"),trait_order)
summary <- fread(file.path(root,"results/gwas_summary.tsv"))[variant=="primary"]
th <- fread(file.path(root,"qa/testing_thresholds.tsv"))
# Independent eligibility check directly from PLINK genotype counts, without
# sourcing the analysis helper that generated the published summary.
counts <- fread(file.path(root,"shared/cohort_counts.frqx"))
setnames(counts,c("chr","SNP","A1","A2","hom1","het","hom2","hap1","hap2","missing"))
stopifnot(!anyDuplicated(counts$SNP),all(counts$hom1+counts$het+counts$hom2+counts$missing==176),
          all(counts$hap1+counts$hap2==0))
counts[,called:=hom1+het+hom2]
counts[,af:=(hom1+het/2)/called]
counts[,`:=`(maf=pmin(af,1-af),call_rate=called/176)]
eligible <- counts[maf>=.10-1e-12 & call_rate>=.90-1e-12,SNP]
M <- length(eligible)
family_cut <- .05/(6*M)
stopifnot(M==1428812L,M==th$eligible_variants,
          abs(family_cut-th$bonferroni_family)<1e-18,all(summary$n_lines==176))
chrom <- data.table(arm=c("2L","2R","3L","3R","4","X"),
                    width=c(23513712,25286936,28110227,32079331,1348131,23542271))
chrom[,offset:=c(0,head(cumsum(width+1000000),-1))]
mp <- list(); qq <- list();checks <- list()
for (tr in trait_order) {
  message("Verifying and rendering primary scan: ",tr)
  scan <- fread(file.path(root,"runs",paste0(tr,"__primary.mlma")),select=c("SNP","bp","p"))
  stopifnot(!anyDuplicated(scan$SNP))
  scan <- scan[match(eligible,SNP)]
  stopifnot(!anyNA(scan$SNP),all(is.finite(scan$p) & scan$p>0 & scan$p<=1))
  # Restore scan order for deterministic downsampling; order does not affect QC.
  setorder(scan,SNP)
  sm <- summary[trait==tr]
  inflation <- median(qchisq(scan$p,df=1,lower.tail=FALSE))/qchisq(.5,df=1)
  stopifnot(nrow(scan)==sm$valid_tests,min(scan$p)==sm$minimum_p,
            sum(scan$p<1e-5)==sm$nominal_hits,abs(inflation-sm$lambda_gc)<1e-10,
            !any(scan$p<.05/M),!any(scan$p<family_cut))
  show <- scan[seq_len(.N)%%60==0 | p<=1e-4]
  stopifnot(all(scan[p<=1e-4,SNP] %in% show$SNP))
  show[,`:=`(trait=tr,arm=sub("_.*$","",SNP),neglogp=-log10(p))]
  show[,position:=(bp+chrom$offset[match(arm,chrom$arm)])/1e6]
  mp[[tr]] <- show
  sorted <- sort(scan$p); n <- length(sorted)
  ranks <- unique(as.integer(c(1:200,round(exp(seq(log(201),log(n),length.out=1000))))))
  qq[[tr]] <- data.table(trait=tr,rank=ranks,n=n,expected=-log10((ranks-.5)/n),
                         observed=-log10(sorted[ranks]),
                         low=-log10(qbeta(.975,ranks,n-ranks+1)),
                         high=-log10(qbeta(.025,ranks,n-ranks+1)))
  checks[[tr]] <- data.table(trait=tr,n_lines=176,eligible_variants=M,valid_tests=n,
    minimum_p=min(scan$p),nominal_hits=sum(scan$p<1e-5),lambda_gc=inflation,
    family_significant=sum(scan$p<family_cut),per_trait_significant=sum(scan$p<.05/M),
    plotted_manhattan=nrow(show),all_p_le_1e4_retained=TRUE,qq_ranks=length(ranks))
  rm(scan,sorted);gc(FALSE)
}
mp <- rbindlist(mp);qq <- rbindlist(qq)
mp[,arm:=factor(arm,levels=chrom$arm)]
mp[,shade:=as.integer(arm)%%2]
write_source(mp,"Figure5_Manhattan_plotted_variants.tsv")
write_source(qq,"Figure5_QQ_order_statistics.tsv")
write_source(rbindlist(checks),"Figure5_independent_scan_checks.tsv")
write_source(th,"Figure5_testing_thresholds.tsv")
write_source(chrom,"Figure5_chromosome_offsets.tsv")

rows <- list()
for (i in seq_along(trait_order)) {
  tr <- trait_order[i]
  pm <- ggplot(mp[trait==tr],aes(position,neglogp)) +
    geom_point(aes(colour=factor(shade)),size=.40,alpha=.70,stroke=0) +
    geom_hline(yintercept=5,linetype="dotted",linewidth=.40,colour="#555555") +
    geom_hline(yintercept=-log10(family_cut),linetype="dashed",linewidth=.45,colour="#262626") +
    scale_colour_manual(values=c("#A0A0A0",blue)) +
    scale_x_continuous(breaks=(chrom$offset+chrom$width/2)/1e6,labels=chrom$arm,
                       expand=expansion(mult=c(.008,.008))) +
    scale_y_continuous(limits=c(0,8.8),breaks=c(0,2,4,6,8),expand=expansion(mult=c(0,0))) +
    labs(title=paste0(LETTERS[i],"  ",labels[tr]),y=expression(-log[10](P)),
         x=if(i==6) "Chromosome arm" else NULL)
  pq <- ggplot(qq[trait==tr],aes(expected,observed)) +
    geom_ribbon(aes(ymin=low,ymax=high),fill="#E5E5E5") +
    geom_abline(slope=1,intercept=0,linetype="dashed",linewidth=.40,colour="#555555") +
    geom_point(colour=blue,size=.48,stroke=0) +
    scale_x_continuous(limits=c(0,7.8),breaks=c(0,2,4,6),expand=expansion(mult=c(0,0))) +
    scale_y_continuous(limits=c(0,7.8),breaks=c(0,2,4,6),expand=expansion(mult=c(0,0))) +
    coord_fixed() + labs(title=sprintf("lambda = %.3f",summary[trait==tr,lambda_gc]),
      x=if(i==6)expression(Expected~-log[10](P)) else NULL,y=expression(Observed~-log[10](P)))
  rows[[i]] <- pm + pq + plot_layout(widths=c(3.7,1.25))
}
p5 <- wrap_plots(rows,ncol=1) + plot_annotation(
  title="Figure 5. Genome-wide association of six primary phenotypes",
  subtitle="176 DGRP lines; call rate >= 90%; MAF >= 10%; 1,428,812 variants per trait",
  caption=paste0("Dotted: exploratory P = 10^-5; dashed: six-trait Bonferroni P = 5.83 x 10^-9.\n",
    "No association passes per-trait or six-trait Bonferroni.\n",
    "QQ shading: pointwise 95% IID-uniform reference; not adjusted for LD."),
  theme=theme(plot.title=element_text(size=12,face="bold"),plot.subtitle=element_text(size=10),
              plot.caption=element_text(size=8.7,hjust=0),plot.margin=margin(7,8,7,7)))
save_figure(p5,"Figure5_primary_GWAS",7.2,10.4,"main")

# Fixed candidates, frozen before this sensitivity analysis: range across six
# alternative fly-QC cohorts, not a confidence interval and not replication.
comp <- fread("manuscript_readiness/candidate_qc/candidate_comparisons.tsv")
rob <- fread("manuscript_readiness/candidate_qc/candidate_robustness_summary.tsv")
stopifnot(nrow(comp)==360,nrow(rob)==60,all(comp$same_direction),all(rob$all_same_direction),
          all(comp[, .N,by=.(trait,SNP)]$N==6),all(comp$alternative_family_significant==FALSE))
independent <- comp[,.(minimum_p=min(p),maximum_p=max(p),schemes_nominal=sum(p<1e-5),
   all_same_direction=all(sign(b)==sign(recorded_beta))),by=.(trait,SNP)]
z <- merge(rob,independent,by=c("trait","SNP"),suffixes=c("","_check"))
stopifnot(all(z$minimum_p==z$minimum_p_check),all(z$maximum_p==z$maximum_p_check),
          all(z$schemes_nominal==z$schemes_nominal_check))
rob[,`:=`(primary_neglogp=-log10(original_p),alternative_low=-log10(maximum_p),
           alternative_high=-log10(minimum_p))]
rob[,trait:=factor(trait,levels=trait_order)]
retained <- rob[,.(original_nominal=sum(recorded_nominal),
   retained_all_six=sum(recorded_nominal & schemes_nominal==6)),by=trait]
stopifnot(sum(retained$original_nominal)==31,sum(retained$retained_all_six)==23)
retained[,label:=sprintf("%d/%d screen hits retained\nin all six alternatives",retained_all_six,original_nominal)]
write_source(rob,"FigureS5_fixed_lead_QC_ranges.tsv")
write_source(comp,"FigureS5_fixed_lead_QC_comparisons.tsv")
write_source(retained,"FigureS5_nominal_retention.tsv")
ps5 <- ggplot(rob,aes(primary_neglogp)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",linewidth=.4,colour="#A0A0A0") +
  geom_hline(yintercept=5,linetype="dotted",linewidth=.4,colour="#555555") +
  geom_vline(xintercept=5,linetype="dotted",linewidth=.4,colour="#555555") +
  geom_linerange(aes(ymin=alternative_low,ymax=alternative_high),colour=blue,linewidth=.7) +
  geom_point(aes(y=primary_neglogp,shape=recorded_nominal),size=1.6,colour=blue,fill="white") +
  geom_text(data=retained,aes(x=3.62,y=6.9,label=label),inherit.aes=FALSE,
            hjust=0,vjust=1,size=3,colour="#333333") +
  scale_shape_manual(values=c(`FALSE`=1,`TRUE`=16)) +
  scale_x_continuous(limits=c(3.5,7),breaks=c(4,5,6,7),expand=expansion(mult=c(0,0))) +
  scale_y_continuous(limits=c(3.5,7),breaks=c(4,5,6,7),expand=expansion(mult=c(0,0))) +
  facet_wrap(~trait,ncol=2,labeller=as_labeller(labels)) +coord_fixed()+
  labs(title="Figure S5. Fixed candidate sensitivity to fly-QC definitions",
    subtitle="10 separated leads per trait; six alternative QC cohorts\nAll effects retain their direction",
    x=expression(Primary~-log[10](P)),y=expression(Alternative~-log[10](P)),
    caption=paste0("Vertical segments: P-value range across six alternatives (not confidence intervals).\n",
       "Points: primary values; filled if P < 10^-5. Dotted: P = 10^-5; dashed: equality.\n",
       "23/31 primary screen hits remain P < 10^-5 in all alternatives.\n",
       "This is sensitivity analysis, not independent replication.")) +
  theme(plot.title=element_text(size=12),plot.caption=element_text(hjust=0,size=9),
        panel.spacing=grid::unit(1.1,"lines"))
save_figure(ps5,"FigureS5_fixed_candidate_QC_sensitivity",7.2,9.0,"supplementary")
message("Saved Figure 5 and Figure S5; all 6 full primary scans independently reconciled.")
