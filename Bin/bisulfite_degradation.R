#!/usr/bin/env RScript
# Title: bisulfite_degradation.R
# Authors: Alexandre Pellan Cheng
# Brief description: Measures degradation after bisulfite treatment (in base pairs)

# Initialize -------------------------------------------------------------------------------------------------
rm(list=ls())

# Load libraries ----------------------------------------------------------------------------------------------------------
library(data.table)
library(ggplot2)
library(parallel)

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
file.list<-list.files(path="../V1/Lengths/")

lengths<-rbindlist(mclapply(X = file.list, FUN = length_stats, mc.cores=51, lengths.path))

ids<-read.table(paste0(lists, "BS_and_nonBS_samples.txt"), header=TRUE)

colnames(ids)[1]<-"sample"

length.df<-merge(ids, lengths, by="sample")
colnames(length.df)<-c("meth_id", "sample", "meth_mean", "meth_median", "meth_sd")

length.df<-merge(length.df, lengths, by="sample")
length.df$decrease<-length.df$length_mean-length.df$meth_mean
length.df$median_decrease<-length.df$length_median-length.df$meth_median
m<-melt(length.df[,c("meth_id", "meth_mean", "length_mean")], id.vars="meth_id")


ggplot(data=m)+geom_boxplot(aes(x=variable, y=value))+geom_point(aes(x=variable, y=value))+geom_line(aes(x=variable, y=value, group=meth_id), color="gray")
