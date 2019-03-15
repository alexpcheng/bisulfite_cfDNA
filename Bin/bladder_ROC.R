#!/usr/bin/env RScript
# Title: Tissue composition by pathology
# Authors: Alexandre Pellan Cheng
# Brief description: Creates figure XYZ

# Initialize -------------------------------------------------------------------------------------------------
rm(list=ls())

# Load libraries ---------------------------------------------------------------------------------------------
library(ggplot2)
library(data.table)
library(pROC)
source("ggplot_n_samples.R")
# Paths ------------------------------------------------------------------------------------------------------
tissue_origin.path <- "../V2/tissues_of_origin/"

# Load files -------------------------------------------------------------------------------------------------
mp.df<-data.frame(fread(paste0(tissue_origin.path, "mp.txt")))
mp.df$group<-gsub("ETP", "Early", mp.df$group)
mp.df$group<-gsub("HLY", "no-UTI", mp.df$group)
mp.df$group<-factor(mp.df$group, levels=c("BKVN", "BKV+/N-", "Normal", "no-UTI", "UTI", "Early"))

# Load total cfDNA
cfDNA<-fread("../lists/sample_cfDNA_extracts")
colnames(cfDNA)[1]<-c("sample")
microbial_fraction<-fread("../V1/read_statistics/BS_treated_samples.txt")
microbial_fraction$human_fraction<-1-(microbial_fraction$microbe_mapped/microbial_fraction$total_reads)
microbial_fraction<-microbial_fraction[, c("sample", "human_fraction")]
cfDNA<-merge(cfDNA, microbial_fraction, by="sample")

cfDNA$humanDNA<-cfDNA$cfDNA_ng_ul*cfDNA$elution_volume/cfDNA$urine_volume_ul*cfDNA$human_fraction

mp.df<-merge(mp.df, cfDNA[, c("sample", "humanDNA"), ], by="sample")
UTI_tested<-droplevels(mp.df[grepl("UTI", mp.df$group), ])

wbc_counts<-fread("../lists/blood_cell_clinical_counts")
UTI_tested<-merge(UTI_tested, wbc_counts, by="sample")
UTI_tested$WBC_count<-as.numeric(UTI_tested$WBC_count)
UTI_tested$RBC_count<-as.numeric(UTI_tested$RBC_count)

cor.test(UTI_tested$humanDNA*UTI_tested$bladder, UTI_tested$RBC_count)

UTI_tested<-UTI_tested[, c("sample", "group", "WBC", "bladder", "kidsum", "humanDNA", "WBC_count", "RBC_count")]
UTI_tested$pyuria<-"No"
UTI_tested$pyuria[UTI_tested$WBC_count>10]<-"Yes"
UTI_tested$colorpyuria<-"blue"
UTI_tested$colorpyuria[UTI_tested$pyuria=="Yes"]<-"red"
UTI_tested$hematuria<-"No"
UTI_tested$hematuria[UTI_tested$RBC_count>4]<-"Yes"

ROC_without_cutoffs<-roc(UTI_tested$group, UTI_tested$bladder*UTI_tested$humanDNA)
ROC_without_cutoffs.df<-data.frame(cbind((1-ROC_without_cutoffs$specificities), ROC_without_cutoffs$sensitivities))
ROC_without_cutoffs.df$state<-"No cutoff"
UTI_only<-UTI_tested[UTI_tested$group=="UTI", ]

pyuria<-UTI_tested[(UTI_tested$pyuria=="Yes" & UTI_tested$group=="UTI") | (UTI_tested$group=="no-UTI"), ]
ROC_pyuria<-roc(pyuria$group, pyuria$bladder*pyuria$humanDNA)
ROC_pyuria.df<-data.frame(cbind((1-ROC_pyuria$specificities), ROC_pyuria$sensitivities))
ROC_pyuria.df$state<-"Pyuria"
ROC_finding_UTI_pyuria<-roc(UTI_only$pyuria, UTI_only$bladder*UTI_only$humanDNA)
ROC_finding_UTI_pyuria.df<-data.frame(cbind((1-ROC_finding_UTI_pyuria$specificities), ROC_finding_UTI_pyuria$sensitivities))
ROC_finding_UTI_pyuria.df$state<-"UTI"

UTI_only<-UTI_tested[UTI_tested$group=="UTI", ]

library(ggbeeswarm)
hematuria<-UTI_tested[(UTI_tested$hematuria=="Yes" & UTI_tested$group=="UTI") | (UTI_tested$group=="no-UTI"), ]

ROC_hematuria<-roc(hematuria$group, hematuria$bladder*hematuria$humanDNA)
ROC_hematuria.df<-data.frame(cbind((1-ROC_hematuria$specificities), ROC_hematuria$sensitivities))
ROC_hematuria.df$state<-"hematuria"
ROC_finding_UTI_hematuria<-roc(UTI_only$hematuria, UTI_only$bladder*UTI_only$humanDNA)
ROC_finding_UTI_hematuria.df<-data.frame(cbind((1-ROC_finding_UTI_hematuria$specificities), ROC_finding_UTI_hematuria$sensitivities))
ROC_finding_UTI_hematuria.df$state<-"UTI"
summary(ROC_finding_UTI_hematuria)

print(ROC_finding_UTI_hematuria$auc)

pdf(file="../Figures/hematuria_bar.pdf",
   width=1.65, height=1.65, paper="special", bg="white", #used to be 1.75 x 1.78
  fonts="Helvetica", colormodel = "cmyk", pointsize=6)
ggplot(UTI_tested, aes(x=group, y=bladder*humanDNA*1000))+
  stat_summary(geom="bar", fun.y=mean, position="dodge", fill="gray", width=0.5)+
  geom_quasirandom(size=0.5, aes(color=hematuria))+scale_color_manual(values=c("blue", "red"))+
  stat_summary(geom="errorbar", fun.data=mean_se, position="dodge", width=0.2)+
  xlab(" ")+
  ylab("Bladder cfDNA (ng/ml)")+
  theme_bw()+
  theme(plot.background=element_blank())+
  theme(legend.position="none")+
  theme(axis.title=element_text(family="Helvetica", size=8),
        axis.text=element_text(family="Helvetica", size=6),
        plot.title=element_blank(), plot.margin = margin(t=5,r=2.5,l=5,b=5))
dev.off()
hematuria.comparison<-rbind(ROC_without_cutoffs.df, ROC_finding_UTI_hematuria.df)

pdf(file="../Figures/hematuria_inset.pdf",
   width=0.7, height=0.7, paper="special", bg="white", #used to be 1.75 x 1.78
  fonts="Helvetica", colormodel = "cmyk", pointsize=6)
ggplot(data=hematuria.comparison, aes(x=X1, y=X2, color=state))+
  geom_segment(aes(x=0, y=0, xend=1, yend=1), color="gray", size=0.5)+
  geom_path(size=0.5)+scale_linetype_discrete()+scale_color_manual(values=c("#1B9E77", "#7570B3"))+
  coord_cartesian(xlim=c(0,1), ylim=c(0,1))+
  scale_x_continuous(breaks=c(0,0.5,1), labels = c(0, 0.5, 1))+
  scale_y_continuous(breaks=c(0,0.5,1), labels=c(0, 0.5, 1))+
  xlab("FPR")+
  ylab("TPR")+
  theme_bw()+
  theme(plot.background=element_blank())+
  theme(legend.position="none")+
  theme(axis.title=element_blank(),
        axis.text.x=element_text(family="Helvetica", size=6, margin=margin(t=0)),
        axis.text.y=element_text(family="Helvetica", size=6, margin=margin(l=-2)),
        axis.ticks.length = unit(1, "pt"),
        plot.margin =margin(t=1, r=1, b=1, l=2),
        plot.title=element_blank(),
        panel.grid = element_blank())
dev.off()