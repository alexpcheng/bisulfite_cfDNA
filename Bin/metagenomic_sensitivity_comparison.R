#!/usr/bin/env RScript
# Title: Metagenomic sample comparison
# Authors: Alexandre Pellan Cheng
# Brief description: Compares GRAMMy outputs
# Load libraries -------------------------------------------------------------------------------------------------------------------------
library(data.table)
library(ggplot2)

# Paths ----------------------------------------------------------------------------------------------------------------------------------
nonbs_samples.path<-"../V1/infections/std_seq/"
bs_samples.path<-"../V1/infections/"
lists<-"../lists/"

# Functions ------------------------------------------------------------------------------------------------------------------------------
source("fancy_scientific.R")

#Get RGE for each detected species
get_infection<-function(filename, nonbs_samples.path, bs_samples.path){
  if (!grepl("MET", filename)){
    file<-data.frame(fread(paste0(nonbs_samples.path, filename, "_BKV.grammy.tab"), header=TRUE))
    f=strsplit(filename, "_")[[1]][1]  
    species.abundance<-file[,c("superkingdom", "species", "RelCoverage")]
    colnames(species.abundance)[3]<-"RelCoverage_std_seq"
    species.abundance<-aggregate(.~species+superkingdom, data=species.abundance, sum)
    species.abundance$std_seq_id<-f
  }
  else{
    file<-data.frame(fread(paste0(bs_samples.path, filename, ".grammy.tab"), header=TRUE))
    f=filename
    species.abundance<-file[,c("superkingdom", "species", "RelCoverage")]
    colnames(species.abundance)[3]<-"RelCoverage_meth_seq"
    species.abundance<-aggregate(.~species+superkingdom, data=species.abundance, sum)
    species.abundance$meth_id<-f
  }
  species.abundance<-species.abundance[(species.abundance$superkingdom!=2759), ]
  species.abundance<-species.abundance[(species.abundance$superkingdom!=2157), ]
  return(species.abundance)
}

# Load samples -------------------------------------------------------------------------------------------------------------------
sams<-fread(paste0(lists, "included_samples.txt"), header=F)
study_ids<-fread(paste0(lists, "BS_and_nonBS_samples.txt"), header=TRUE)

wgbs.rge<-rbindlist(lapply(X = as.list(study_ids$meth_id), FUN = get_infection, nonbs_samples.path, bs_samples.path))
std.rge<-rbindlist(lapply(X=as.list(study_ids$std_seq_id), FUN=get_infection, nonbs_samples.path, bs_samples.path))

rge<-merge(wgbs.rge, study_ids, by="meth_id")
rge<-merge(rge, std.rge, by=c("std_seq_id", "species", "superkingdom"), all=TRUE)
rge<-merge(rge[,-c("meth_id")], study_ids, by="std_seq_id")

rge$missing<-apply(X = rge, MARGIN = 1, FUN = function(x) if (any(is.na(x))) return(0) else return(1))
cor.test(rge$RelCoverage_meth_seq, rge$RelCoverage_std_seq, method="spearman")
save_eps<-TRUE

rge$diff<-abs(rge$RelCoverage_meth_seq/rge$RelCoverage_std_seq)

if(save_eps){pdf(file="../Figures/RGE_std_vs_meth.pdf",
                 width=1.625, height=1.625, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6)}
ggplot(data=rge)+
  geom_point(aes(x=RelCoverage_std_seq, y=RelCoverage_meth_seq, color=factor(superkingdom)), size=0.5, fill="NA", pch=21)+
  scale_x_log10(labels=fancy_scientific)+
  scale_y_log10(labels=fancy_scientific)+
  coord_cartesian(xlim=c(10^-3:15^5), ylim=c(10^-3:15^5))+
  scale_color_manual(values=c("#00B247", "purple"))+
  scale_shape_manual(values=c(4,1))+
  xlab(" ")+
  ylab(" ")+
  theme_bw()+
  theme(plot.background=element_blank())+
  theme(legend.position="none")+
  theme(axis.title.x=element_text(family="Helvetica", size=6),
        axis.title.y=element_text(family="Helvetica", size=6),
        axis.text=element_text(family="Helvetica", size=6),
        plot.title=element_text(family="Helvetica", size=6))
if(save_eps){dev.off()}

