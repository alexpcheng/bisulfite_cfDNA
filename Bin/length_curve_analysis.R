#!/usr/bin/env RScript
# Title: [ENTER TITLE]
# Authors: [ENTER AUTHORS]
# Brief description: [ENTER DESCRIPTION]

# Initialize -------------------------------------------------------------------------------------------------
rm(list=ls())

# Load libraries ----------------------------------------------------------------------------------------------------------
library(data.table)
library(ggplot2)
library(parallel)
source("ggplot_n_samples.R")

# Paths -------------------------------------------------------------------------------------------------------------------
lists<- "../lists/"
lengths.path<-"../V1/Lengths/"

# Functions ---------------------------------------------------------------------------------------------------------------
length_stats<-function(filename, lenghts.path){
  sample<-strsplit(filename, "[.]")[[1]][1]
  sam<-fread(paste0(lengths.path,filename))
  length_mean<-mean(sam$V1)
  length_median<-median(sam$V1)
  std_dev<-sd(sam$V1)
  df<-data.frame(sample, length_mean, length_median, std_dev)
  return(df)
}

#Load data ----------------------------------------------------------------------------------------------------------------
sample_groups<-data.frame(fread(paste0(lists, "sample_groups.txt"), header=FALSE))
colnames(sample_groups)<-c("sample", "group")
file.list<-list.files(path =lengths.path, pattern = "lengths")

lengths<-rbindlist(mclapply(X = file.list, FUN = length_stats, mc.cores=5, lengths.path))

lengths<-lengths[grepl("MET", lengths$sample), ]

mp.df<-fread("../V2/tissues_of_origin/mp.txt")
cfDNA<-fread("../lists/sample_cfDNA_extracts")
microbial_reads<-fread("../V1/microbial_read_fraction.txt")
cfDNA<-merge(cfDNA, microbial_reads, by="meth_id")

cfDNA$humanDNA<-cfDNA$elution_volume*cfDNA$cfDNA_ng_ul/cfDNA$urine_volume_ul*cfDNA$human_fraction
colnames(cfDNA)[1]<-"sample"
lengths<-merge(lengths, mp.df[, c("sample", "WBC", "kidsum")], by="sample")
lengths<-merge(lengths, cfDNA[, c("sample", "humanDNA")], by="sample")
lengths<-merge(lengths, sample_groups, by="sample")
lengths$group<-gsub("HLY", "No UTI", lengths$group)
lengths$group<-gsub("vir-/inf-", "Normal", lengths$group)
lengths$group<-gsub("ETP", "Early", lengths$group)
lengths$group<-gsub("vir[+]/inf[+]|vir[+]/inf-", "BKV+/N-", lengths$group)

lengths$group<-factor(lengths$group, levels=c("BKVN", "BKV+/N-", "Normal", "No UTI", "UTI", "Early"))

cor.test(lengths$length_mean, lengths$WBC, method="spearman")
wilcox.test(lengths[grepl("UTI|No UTI", lengths$group), ]$length_mean~lengths[grepl("UTI|No UTI", lengths$group), ]$group)

mean(lengths[lengths$group=="UTI", ]$length_mean)
mean(lengths[lengths$group=="No UTI", ]$length_mean)

ggplot(data=lengths)+geom_boxplot(aes(x=group, y=length_mean))+geom_point(aes(x=group, y=length_mean))+
  theme_bw()+
  xlab("Group")+ylab("Mean length (bp)")


save_eps<-T
if(save_eps){pdf(file="../Figures/lengths.pdf",
                 width=0.65, height=0.63, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6)}
ggplot(data=lengths[grepl("UTI", lengths$group), ], aes(x=group, y=length_mean, fill=group))+
  geom_boxplot(outlier.size = 0.1, color="black", size=0.25)+geom_point(position=position_dodge(width=0.75), size=0.1)+
  geom_vline(xintercept=4.5, color="black", linetype="dashed", size=0.25)+
  scale_fill_manual(values=c("#984EA3", "#FF7F00"))+
  scale_y_continuous(breaks=c(90, 110))+
  ylab("Mean fragment\nlength (bp)")+
  theme_bw()+
  theme(plot.background=element_blank())+
  theme(legend.position="none")+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_text(family="Helvetica", size=6, margin = margin(t=-1)),
        axis.text=element_text(family="Helvetica", size=6),
        axis.text.x=element_blank(),
        axis.ticks.length = unit(1, "pt"),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(margin=margin(l=-1)),
        plot.margin =margin(t=1, r=1, b=1, l=1),
        plot.title=element_blank(),
        panel.grid = element_blank())
if(save_eps){dev.off()}

if(save_eps){pdf(file="../Figures/length_vs_WBC.pdf",
                 width=1.65, height=1.65, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6)}
ggplot(data=lengths, aes(x=WBC, y=length_mean))+
  geom_point(size=0.75)+geom_smooth(method="lm", color="#696969", se=FALSE)+
  scale_fill_manual(values=c("blue", "red"))+
  ylab("Mean fragment length (bp)")+ylim(c(50, 125))+
  xlab("Leukocyte proportion")+
  theme_bw()+
  theme(plot.background=element_blank())+
  theme(legend.position="none")+
  theme(axis.title=element_text(family="Helvetica", size=8),
        axis.text=element_text(family="Helvetica", size=6),
        plot.margin =margin(t=5, r=1, b=5, l=5),
        plot.title=element_blank())
if(save_eps){dev.off()}

if(save_eps){pdf(file="../Figures/lengths_all.pdf",
                 width=3, height=3, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6)}
ggplot(data=lengths, aes(x=group, y=length_mean))+
  geom_boxplot(outlier.size = 1, color="black")+geom_point(position=position_dodge(width=0.75), size=1)+
  n_samples(m, x=lengths$group, y=rep(0, nrow(lengths)), y_nudge=83, fontsize=6)+
  geom_vline(xintercept=5.5, color="black", linetype="dashed", size=0.25)+
  ylab("Mean fragment\nlength (bp)")+xlab("Group")+
  theme_bw()+
  theme(plot.background=element_blank())+
  theme(legend.position="none")+
  theme(axis.title=element_text(family="Helvetica", size=8),
        axis.text=element_text(family="Helvetica", size=6),
        panel.grid = element_blank())
if(save_eps){dev.off()}

