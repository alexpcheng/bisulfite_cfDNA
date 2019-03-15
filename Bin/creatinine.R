#!/usr/bin/env RScript
# Title: creatinine.R
# Authors: Alexandre Pellan Cheng
# Brief description: Looks at correlation between creatinine and kidney cfDNA

# Initialize -------------------------------------------------------------------------------------------------
rm(list=ls())

# Load libraries ---------------------------------------------------------------------------------------------
library(ggplot2)
library(data.table)
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

creat<-fread("../lists/creatinine.list", header=TRUE)
creat$creatinine<-as.numeric(creat$creatinine)

mp.df<-merge(mp.df, creat, by="sample")

pdf(file="../Figures/creatinine_v3.pdf",
    width=1.65, height=1.65, paper="special", bg="white",
    fonts="Helvetica", colormodel = "cmyk", pointsize=6)
ggplot(data=mp.df)+geom_point(aes(y=kidsum*humanDNA*1000, x=creatinine, color=group, shape=group), size=0.5)+
  scale_color_manual(values= c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628"))+
  scale_shape_manual(values= c(19, 17, 3, 5, 4, 8))+
  ylab("Kidney cfDNA (ng/ml)")+ylim(values=c(0,65))+
  xlab("Creatinine (mg/dL)")+
  theme_bw()+coord_cartesian()+
  theme(plot.background=element_blank())+
  theme(legend.position="none")+
  theme(axis.title=element_text(family="Helvetica", size=8),
        axis.text=element_text(family="Helvetica", size=6))
dev.off()