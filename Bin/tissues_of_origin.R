#!/usr/bin/env RScript
# Title: Tissue of Origin Measurements
# Authors: Alexandre Pellan Cheng
# Brief description: Creates a text file with the tissue of origin measurements (aka mixing parameters (mp)) for
# sample. References are loaded and have already assumed to be prepped (via references.500.rds)

# Initialize -------------------------------------------------------------------------------------------------
rm(list=ls())
# Load libraries ---------------------------------------------------------------------------------------------
library(Matrix)
library(matrixcalc)
library(stats)
library(data.table)
library(limSolve)
library(parallel)

# Source custom function ------------------------------------------------------------------------------------- 
source("tissues_of_origin_functions.R")

# Paths ------------------------------------------------------------------------------------------------------
samples.path <- "../V1/binned_samples/"
references.path <-"../Methylation_References_snaked/MethylMatrix/" #usually "../Methylation_References/
lists.path <- "../lists/"
donor_fractions.path <- "../V2/donor_fractions/"
tissue_origin.path <- "../V2/tissues_of_origin/"

dir.create(tissue_origin.path, showWarnings = F)
# Load samples and references ---------------------------------------------------------------------------------
file_list = list.files(path = samples.path, pattern="*MET")

# Function to read and slightly modify samples
process_samples<-function(filename, samples.path){
  tiled<-fread(paste0(samples.path, filename))
  colnames(tiled)<-c("chr", "start", "end", "pmeth")
  tiled<-tiled[complete.cases(tiled),]
  tiled$pmeth<-tiled$pmeth/100
  colnames(tiled)[4]<-strsplit(filename, ".filtered")[[1]][1]
  return(tiled)
}

permeth_samples<-mclapply(file_list, FUN=process_samples, samples.path, mc.cores=10)
#reference_tissues<-fread(paste0(references.path, "MM.binned"))
reference_tissues<-fread(paste0(references.path, "MethylMatrix_binned"))
#coln<-fread("../Methylation_References/groups_and_reference_key.list", header=F)
#colnames(reference_tissues)<-c("chr", "start", "end", coln$V2)

coln<-fread("../Methylation_References_snaked/lookup_table.txt", header=F)
colnames(reference_tissues)<-c("chr", "start", "end", coln$V3)
# Initialize variables ----------------------------------------------------------------------------------------
num.tissues<-ncol(reference_tissues)-3 # 3 columns are for chromosome, start and end
num.samples<-length(file_list)

sample_groups<-fread(paste0(lists.path, "sample_groups.txt"), header = F)
colnames(sample_groups)<-c("sample", "group")
sample_groups$group<-gsub("vir-/inf-", "Normal", sample_groups$group)
sample_groups$group<-gsub("vir[+]/inf-", "BKV+/N-", sample_groups$group)
sample_groups$group<-gsub("vir[+]/inf[+]", "BKV+/N-", sample_groups$group)
donor_fractions<-fread(paste0(donor_fractions.path, "donor_fractions.txt"), header=T)
colnames(donor_fractions)<-c("sample", "methDF")


# Table preparation ------------------------------------------------------------------------------------------
merr<-function(sam, ref){
  z=merge(ref, sam, by=c("chr", "start", "end"))
  z=z[complete.cases(z)]
  return(z)
}

reference_and_samples<-mclapply(FUN=merr,
                                X=permeth_samples, reference_tissues,
                                mc.cores = 10)

# Measuring tissue of orign ----------------------------------------------------------------------------------
# Create a list where each element contains the tissues of origin and the error measurement

tissues_of_origin.list<-mclapply(reference_and_samples, FUN=mp_function,other=TRUE,
                                 nt=num.tissues, sum_to_one=FALSE, mc.cores=10)
tissues_of_origin<-c()
for (i in 1:num.samples){
  tissues_of_origin<-c(tissues_of_origin, tissues_of_origin.list[[i]])
}

mp.tmp<-Reduce(x=tissues_of_origin, f=function(x,y) merge(x,y, by="tissue", all=TRUE))
mp.df<-as.data.frame(t(mp.tmp[,2:ncol(mp.tmp)]))
colnames(mp.df)<-mp.tmp$tissue
mp.df$sample<-colnames(mp.tmp)[2:(num.samples+1)]
mp.df$kidsum<-rowSums(mp.df[, grep("kidney|podocyte|mesangial", colnames(mp.df))])
mp.df$lympho<-rowSums(mp.df[,grepl("NK|TCell|BCell", colnames(mp.df))])
mp.df$myelo<-rowSums(mp.df[,grepl("macrophage|monocyte|dendritic|eosonophil|neutrophil", colnames(mp.df))])
mp.df$WBC<-mp.df$lympho+mp.df$myelo
mp.df$bladder<-rowSums(mp.df[, grep("bladder", colnames(mp.df))])
mp.df<-merge(mp.df, donor_fractions, by="sample", all.x=TRUE)
mp.df<-merge(mp.df, sample_groups, by="sample")
fwrite(mp.df, file = paste0(tissue_origin.path, "mp.txt"), quote=FALSE, sep="\t")

