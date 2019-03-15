#!/usr/bin/env RScript
# Title: UTI_sensitivie and specificity
# Authors: Alexandre Pellan Cheng
# Brief description: Measures accuracy of detection for UTI detection

# Initialize -------------------------------------------------------------------------------------------------
rm(list=ls())

# Load libraries ----------------------------------------------------------------------------------------------------------
library(data.table)
library(ROCR)
library(plyr)
library(ggplot2)

# Paths -------------------------------------------------------------------------------------------------------------------
lists<- "../lists/"
infection.path<-"../V1/infections/"
# Functions ---------------------------------------------------------------------------------------------------------------
get_common_uropath<-function(filename, infection.path){
  file<-data.frame(fread(paste0(infection.path, filename, ".grammy.tab"), header=TRUE))
  species.abundance<-file[,c("species", "RelCoverage")]
  colnames(species.abundance)[2]<-"RelCoverage_meth_seq"
  common_uropath.abundance<-species.abundance[which(species.abundance$species==562 | # ecoli
                                                      species.abundance$species==1351 | # e faecalis
                                                      species.abundance$species==287 | # p aeruginosa
                                                      species.abundance$species==571 | # k oxytoca
                                                      species.abundance$species==1352),] # e faecium
  # staph genus exception
  genus.abundance<-file[,c("genus", "RelCoverage")]
  colnames(genus.abundance)[2]<-"RelCoverage_meth_seq"
  common_uropath.genus<-genus.abundance[which(genus.abundance$genus==1279),]
  common_uropath.genus<-aggregate(.~genus, data=common_uropath.genus, sum)
  common_uropath.genus$meth_id<-filename
  
  common_uropath.abundance<-aggregate(.~species, data=common_uropath.abundance, sum)
  common_uropath.abundance$meth_id<-filename
  common_uropath.abundance<-rbind.fill(common_uropath.abundance, common_uropath.genus)
  
  bugs<-c(562, 1351, 287, 571, 1352)
  for (i in bugs){
    print(i)
    if (length(which(common_uropath.abundance$species==i))==0){
      df<-data.frame(matrix(ncol=4, nrow=1))
      colnames(df)<-colnames(common_uropath.abundance)
      df$meth_id<-filename
      df$species<-i
      df$RelCoverage_meth_seq<-0
      common_uropath.abundance<-rbind(common_uropath.abundance, df)
    }
  }
  if (length(which(common_uropath.abundance$genus==1279))==0){
    df<-data.frame(matrix(ncol=4, nrow=1))
    colnames(df)<-colnames(common_uropath.abundance)
    df$meth_id<-file
    df$species<-1279
    df$RelCoverage_meth_seq<-0
    common_uropath.abundance<-rbind(common_uropath.abundance, df)
  }
  return(common_uropath.abundance)
}

#Load data ----------------------------------------------------------------------------------------------------------------
sams<-fread(paste0(lists, "included_samples.txt"), header=F)
urine_culture<-data.frame(fread(paste0(lists, "urine_culture")))

common_uropath<-rbindlist(lapply(X = as.list(as.data.frame(t(sams))), FUN = get_common_uropath, infection.path))

common_uropath<-merge(common_uropath, urine_culture, by="meth_id")
common_uropath<-common_uropath[common_uropath$urine_culture=="Yes", ]

ecoli_thresh<-min(common_uropath[common_uropath$primary=="Escherichia coli" & common_uropath$species==562, ]$RelCoverage_meth_seq)
paeru_thresh<-min(common_uropath[common_uropath$primary=="Pseudomonas aeruginosa" & common_uropath$species==287, ]$RelCoverage_meth_seq)
efaecalis_thresh<-min(common_uropath[common_uropath$primary=="Enterococcus faecalis" & common_uropath$species==1351, ]$RelCoverage_meth_seq)
efaecium_thresh<-min(common_uropath[common_uropath$primary=="Enterococcus faecium" & common_uropath$species==1352, ]$RelCoverage_meth_seq)
koxy_thresh<-min(common_uropath[common_uropath$primary=="Klebsiella oxytoca" & common_uropath$species==571, ]$RelCoverage_meth_seq)
cons_thresh<-min(common_uropath[common_uropath$secondary=="Staphylococcus (coagulase negative)" & common_uropath$genus==1279, ]$RelCoverage_meth_seq)


detect_pathogen<-function(bug, common_uropath){
  microbe=bug[1]
  taxid=bug[2]
  thresh=as.numeric(bug[3])-0.0000001 #accounts for tol issues
  if (microbe != "Staphylococcus (coagulase negative)"){
    TP<-nrow(common_uropath[common_uropath$primary==microbe & 
                                common_uropath$species==taxid & 
                                common_uropath$RelCoverage_meth_seq>=thresh, ])
    FP<-nrow(common_uropath[common_uropath$primary!=microbe & 
                                common_uropath$species==taxid & 
                                common_uropath$RelCoverage_meth_seq>=thresh, ])
    FN<-nrow(common_uropath[common_uropath$primary==microbe & 
                                common_uropath$species==taxid & 
                                common_uropath$RelCoverage_meth_seq<thresh, ])
    TN<-nrow(common_uropath[common_uropath$primary!=microbe & 
                                common_uropath$species==taxid & 
                                common_uropath$RelCoverage_meth_seq<thresh, ])
    print(TP)
  }
  else{
    TP<-nrow(common_uropath[common_uropath$secondary==microbe & 
                                common_uropath$genus==taxid & 
                                common_uropath$RelCoverage_meth_seq>=thresh, ])
    FP<-nrow(common_uropath[common_uropath$secondary!=microbe & 
                                common_uropath$genus==taxid & 
                                common_uropath$RelCoverage_meth_seq>=thresh, ])
    FN<-nrow(common_uropath[common_uropath$secondry==microbe & 
                                common_uropath$genus==taxid & 
                                common_uropath$RelCoverage_meth_seq<thresh, ])
    TN<-nrow(common_uropath[common_uropath$secondary!=microbe & 
                                common_uropath$genus==taxid & 
                                common_uropath$RelCoverage_meth_seq<thresh, ])
    }
  df<-data.frame(microbe, TP, FP, FN, TN)
  return(df)
}

detect_pathogen<-function(bug, common_uropath){
  microbe<-bug[1]
  taxid<-bug[2]
  thresh<-as.numeric(bug[3])-0.00000001
  
  if(microbe!="Staphylococcus (coagulase negative)"){
    c_u<-common_uropath[common_uropath$species==taxid, ]
    c_u$predict<-NA
    c_u$predict[c_u$RelCoverage_meth_seq>=thresh]<-microbe
    c_u$predict[c_u$RelCoverage_meth_seq<thresh]<-paste0(microbe, "NEG")
    c_u$ref<-paste0(microbe, "NEG")
    c_u$ref[c_u$primary==microbe | c_u$secondary==microbe]<-microbe
  }
  else{
    c_u<-common_uropath[common_uropath$genus==taxid, ]
    c_u$predict<-NA
    c_u$predict[c_u$RelCoverage_meth_seq>=thresh]<-microbe
    c_u$predict[c_u$RelCoverage_meth_seq<thresh]<- paste0(microbe, "NEG")
    c_u$ref<-paste0(microbe, "NEG")
    c_u$ref[c_u$primary==microbe | c_u$secondary==microbe]<-microbe
  }
  return(c_u[, c("predict", "ref")])
}

paeru<-detect_pathogen(c("Pseudomonas aeruginosa", 287, paeru_thresh), common_uropath)
efaecalis<-detect_pathogen(c("Enterococcus faecalis", 1351, efaecalis_thresh), common_uropath)
ecoli<-detect_pathogen(c("Escherichia coli", 562, ecoli_thresh), common_uropath)
koxy<-detect_pathogen(c("Klebsiella oxytoca", 571, koxy_thresh), common_uropath)
efaecium<-detect_pathogen(c("Enterococcus faecium", 1352, efaecium_thresh), common_uropath)
cons<-detect_pathogen(c("Staphylococcus (coagulase negative)", 1279, cons_thresh), common_uropath)

all_values<-rbind(paeru, efaecalis, ecoli, koxy, efaecium, cons)

table(all_values$predict, all_values$ref)
caret::confusionMatrix(table(all_values$predict, all_values$ref))
