#!/usr/bin/env RScript
# Title: BKV_qPCR
# Authors: Alexandre Pellan Cheng
# Brief description: Looks at correlation between plasma BKV abundance and BKV relative genomic abundance through sequencing of urinary cfDNA

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

qpcr<-fread("../lists/bkv_qpcr.list", header=FALSE, col.names=c("samples", "qPCR"))
qpcr$qPCR<-as.numeric(qpcr$qPCR)

BKV_rges<-merge(BKV_rges, qpcr, by="samples")
ggplot(data=BKV_rges)+geom_point(aes(x=qPCR, y=BK_rge, color=group))+scale_x_log10()
ggplot(data=BKV_rges)+geom_violin(aes(x=group, y=BK_rge))+scale_y_log10()
cor.test(BKV_rges$qPCR, BKV_rges$BK_rge, method="spearman")
