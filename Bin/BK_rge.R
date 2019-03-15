#!/usr/bin/env RScript
# Title: BK relative genomic abundance
# Authors: Alexandre Pellan Cheng
# Brief description: Creates violin plot of BK relative genomic abundance

# Initialize -------------------------------------------------------------------------------------------------
rm(list=ls())

# Load libraries ---------------------------------------------------------------------------------------------
library(ggplot2)
library(data.table)
library(scales)

# Paths ------------------------------------------------------------------------------------------------------
infections.path <- "../V1/infections/"
samples.list<- "../lists/included_samples.txt"
samples.groups<- "../lists/sample_groups.txt"

# Functions --------------------------------------------------------------------------------------------------
source("fancy_scientific.R")
source("ggplot_n_samples.R")
get_BK<-function(filename, infections.path){
  file<-fread(paste0(infections.path, filename, ".grammy.tab"), header=TRUE)
  BK_rge<-file[grepl("BK_polyomavirus", file$Name),]$RelCoverage
  if (length(BK_rge)==0)
    BK_rge<-0
  df<-data.frame(filename, BK_rge)
  return(df)
}

sample_groups<-fread(samples.groups, header=F)
colnames(sample_groups)<-c("samples", "group")
sams<-fread(samples.list, header=F)
BKV_rges<-rbindlist(apply(X = sams, MARGIN=1, FUN = get_BK, infections.path))
colnames(BKV_rges)<-c("samples", "BK_rge")
BKV_rges<-merge(BKV_rges, sample_groups, by="samples")
BKV_rges$group<-gsub("vir-/inf-", "Normal", BKV_rges$group)
BKV_rges$group<-gsub("vir[+]/inf-", "BKV+/N-", BKV_rges$group)
BKV_rges$group<-gsub("vir[+]/inf[+]", "BKV+/N-", BKV_rges$group)
BKV_rges$group<-gsub("ETP", "Early", BKV_rges$group)
BKV_rges$group<-gsub("HLY", "no-UTI", BKV_rges$group)

BKV_rges$group<-factor(BKV_rges$group, levels=c("BKVN", "BKV+/N-", "Normal", "no-UTI","UTI", "Early"))

save_eps<-T

if(save_eps){pdf(file="../Figures/BKV_rge.pdf",
                 width=1.7, height=1.7, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6)}
ggplot(data=BKV_rges, aes(x=group, y=BK_rge+0.01))+geom_violin(color="black", size=0.5, scale = "width")+
  geom_point(size=0.5)+
  n_samples(BKV_rges, x=BKV_rges$group, y=rep(0.01, nrow(BKV_rges)), y_nudge=-0.4, fontsize=6)+
  scale_y_log10(labels=fancy_scientific)+
  ylab("BKV abundance (RGE)")+
  theme_bw()+
  theme(plot.background=element_blank())+
  theme(legend.position="none")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(family="Helvetica", size=6),
        axis.text=element_text(family="Helvetica", size=6),
        axis.text.x = element_text(family="Helvetica", size=6, angle=40, hjust=1),
        plot.title=element_text(family="Helvetica", size=6))
if(save_eps){dev.off()}
