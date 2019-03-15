#!/usr/bin/env RScript
# Title: UTI_detection.R
# Authors: Alexandre Pellan Cheng
# Brief description: creates figure

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
                                                      #species.abundance$species==573 | # k. pneumoniae
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
    common_uropath.abundance<-rbind(common_uropath.abundance, df)
  }
  return(common_uropath.abundance)
}

# Get data for ROC
ROC<-function(taxid, tested_uropath){
  if (taxid==562){
    coln<-10
  }
  if (taxid==1351){
    coln<-9
  }
  if (taxid==287){
    coln<-8
  }
  if (taxid==573){
    coln<-12
  }
  if (taxid==571){
    coln<-11
  }
  if (taxid==1352){
    coln<-13
  }
  if (taxid==1279){
    coln<-14
  }
  tested_uropath<-data.frame(tested_uropath)
  subset_pathogen<-tested_uropath[(which(tested_uropath$species==taxid | tested_uropath$genus==taxid)),]
  
  pred<-prediction(subset_pathogen$RelCoverage_meth_seq,
                          subset_pathogen[, coln])
  
  perf<-performance(pred, measure= "tpr", x.measure="fpr")
  print(tested_uropath)
  print(perf)
  df<-data.frame(cbind(perf@x.values[[1]], perf@y.values[[1]]))
  colnames(df)<-c("FP", "TP")
  df$taxid<-taxid
  return(df)
}

#Load data ----------------------------------------------------------------------------------------------------------------
sams<-fread(paste0(lists, "included_samples.txt"), header=F)
urine_culture<-data.frame(fread(paste0(lists, "urine_culture")))

common_uropath<-rbindlist(lapply(X = as.list(as.data.frame(t(sams))), FUN = get_common_uropath, infection.path))

common_uropath<-merge(common_uropath, urine_culture, by="meth_id")

test<-common_uropath

test<-test[grepl("Yes", test$urine_culture),]
test$bug<-test$species
test$bug[is.na(test$bug)]<-"Co.NS"
test$bug<-gsub("573", "K. pneumoniae", test$bug)
test$bug<-gsub("571", "K. oxytoca", test$bug)
test$bug<-gsub("562", "E.coli", test$bug)
test$bug<-gsub("287", "P. aeruginosa", test$bug)
test$bug<-gsub("1352", "E. faecium", test$bug)
test$bug<-gsub("1351", "E. faecalis", test$bug)

test$meth_id<-factor(test$meth_id, 
                     levels=c("MET-1-18", "MET-1-30", "MET-1-36", "MET-1-45", 
                              "MET-1-03", "MET-1-07", "MET-1-32", "MET-1-47",
                              "MET-1-34",
                              "MET-1-25", 
                              "MET-1-26", "MET-1-01",
                              "MET-1-10", "MET-1-11", "MET-1-14", "MET-1-16",
                              "MET-1-20", "MET-1-23", "MET-1-24", "MET-1-31", "MET-1-33",
                              "MET-1-37", "MET-1-39", "MET-1-44"))


test$bug<-factor(test$bug, levels=c("Co.NS","P. aeruginosa", "K. oxytoca", "E. faecium", "E. faecalis",  "E.coli"))
test$ceiled<-test$RelCoverage_meth_seq
test$ceiled[test$ceiled>=2]<-2

test$culture<-paste(test$primary, test$secondary, sep="+")
test$culture<-gsub("[+]N/A", "", test$culture)
test$culture<-gsub("IG", "NEG", test$culture)
test$culture<-gsub("Staphylococcus [(]coagulase negative[)]", "Co.NS", test$culture)

test$culture<-gsub("Klebsiella oxytoca", "K. oxytoca", test$culture)
test$culture<-gsub("Escherichia coli", "E. coli", test$culture)
test$culture<-gsub("Pseudomonas aeruginosa", "P. aeruginosa", test$culture)
test$culture<-gsub("Enterococcus faecium", "E. faecium", test$culture)
test$culture<-gsub("Enterococcus faecalis", "E. faecalis", test$culture)

test$culture<-factor(test$culture, levels=c("E. coli", "E. faecalis", "E. faecium", "K. oxytoca", "P. aeruginosa", "E. faecalis+Co.NS", "NEG"))
save_eps<-T

if(save_eps){pdf(file="../Figures/UTI_RGE.pdf",
                 width=3.325, height=1.15, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6)}
ggplot(data=test)+geom_tile(aes(x=meth_id, y=bug, fill=ceiled), color="black")+
  facet_grid(~culture, space="free", scales="free")+theme_bw()+
  theme(strip.text = element_blank(), strip.background = element_blank(), panel.spacing = unit(3, "pt"))+
  scale_fill_gradient(low = "white", high="red")+
  theme(plot.background=element_blank())+
  theme(legend.position="none")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(family="Helvetica", size=6, margin = margin(t=5, r=0, b=5, l=0)),
        plot.title=element_text(family="Helvetica", size=6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(t=5, r=5, b=5, l=0),
        axis.ticks.y=element_blank())
if(save_eps){dev.off()}
