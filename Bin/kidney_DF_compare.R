#!/usr/bin/env RScript
# Title: Kidney vs Donor fraction
# Authors: Alexandre Pellan Cheng
# Brief description: Creates figure XYZ

# Initialize -------------------------------------------------------------------------------------------------
rm(list=ls())
# Load libraries ---------------------------------------------------------------------------------------------
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(data.table)

# Paths ------------------------------------------------------------------------------------------------------
tissue_origin.path <- "../V2/tissues_of_origin/"
lists<- "../lists/"

# Load tissues of origin -------------------------------------------------------------------------------------
mp.df<-fread(paste0(tissue_origin.path, "mp.txt"))

# Load donor fractions from non bisulfite treated sequencing ------------------------------------------------
donor_fractions_std<-fread(paste0(lists, "donor_fractions_std.list"))
colnames(donor_fractions_std)<-c("sample", "std_df")

#Format for plotting
mp.df<-mp.df[!grepl("MET-1-10|MET-1-18", mp.df$sample),] # Remove samples that also received a bone marrow transplant

mp.df<-merge(mp.df[, c("sample", "kidsum", "methDF", "group", "rel_error")], donor_fractions_std, by="sample", all.x=TRUE)
mp.df$methDF<-mp.df$methDF/100
mp.df$std_df<-mp.df$std_df/100
mp.df$group<-factor(mp.df$group, levels=c("BKVN", "BKV+/N-", "ETP", "HLY", "UTI"))
mp.df$df_error<-abs(mp.df$methDF-mp.df$std_df)
non.sexmm<-mp.df[which(is.na(mp.df$methDF)),c("kidsum", "group", "methDF", "rel_error")]
non.sexmm$methDF[which(is.na(non.sexmm$methDF))]<-"N.S."

corelation<-cor.test(mp.df[!is.na(mp.df$methDF),]$kidsum, mp.df[!is.na(mp.df$methDF), ]$methDF, method="spearman")
# Plot data -------------------------------------------------------------------------------------------------

save_eps=TRUE

if(save_eps){pdf(file="../Figures/kidney_total_errorbars.pdf",
                 width=1.825, height=1.825, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats = FALSE)}
ggplot(data=mp.df)+geom_abline(intercept=0)+
  #geom_errorbar(aes(x=kidsum, ymin=(methDF-df_error/2), ymax=(methDF+df_error/2), width=0.03), size=0.2)+
  #geom_errorbarh(aes(xmin=(kidsum-rel_error/2), xmax=(kidsum+rel_error/2), y=methDF, height=0.03), size=0.2)+
  geom_point(aes(y=methDF, x=kidsum), color="black", size=1.0, stroke=0.2)+
  xlab("Kidney fraction")+
  ylab("Donor fraction")+
  theme_bw()+coord_cartesian(ylim=c(0,1), xlim=c(0,1))+
  theme(plot.background=element_blank())+
  theme(legend.position="none")+
  theme(axis.title.x=element_text(family="Helvetica", size=8), axis.title.y=element_text(family="Helvetica", size=8),
        axis.text=element_text(family="Helvetica", size=6),
        plot.title=element_text(family="Helvetica", size=6))
if(save_eps)(dev.off())
