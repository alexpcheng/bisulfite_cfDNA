#!/usr/bin/env RScript
# Title: Tissue composition by pathology
# Authors: Alexandre Pellan Cheng
# Brief description: Creates figure XYZ

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

#Prep WBC
mp.tissue_groups<-mp.df[,c("sample", "group"), drop=FALSE]
mp.tissue_groups$macrophage<-rowSums(mp.df[,grepl("macrophage", colnames(mp.df)), drop=FALSE])
mp.tissue_groups$dendritic<-rowSums(mp.df[,grepl("dendritic", colnames(mp.df)), drop=FALSE])
mp.tissue_groups$eosonophil<-rowSums(mp.df[,grepl("eosonophil", colnames(mp.df)), drop=FALSE])
mp.tissue_groups$neutrophil<-rowSums(mp.df[,grepl("neutrophil", colnames(mp.df)), drop=FALSE])
mp.tissue_groups$NKCell<-rowSums(mp.df[,grepl("NKCell", colnames(mp.df)), drop=FALSE])
mp.tissue_groups$monocyte<-rowSums(mp.df[,grepl("monocyte", colnames(mp.df)), drop=FALSE])
mp.tissue_groups$TCell<-rowSums(mp.df[,grepl("TCell", colnames(mp.df)), drop=FALSE])
mp.tissue_groups$BCell<-rowSums(mp.df[,grepl("BCell", colnames(mp.df)), drop=FALSE])

mp.tissue_groups$bladder<-rowSums(mp.df[,grepl("bladder", colnames(mp.df)), drop=FALSE])
mp.tissue_groups$colon<-rowSums(mp.df[,grepl("colon", colnames(mp.df)), drop=FALSE])
mp.tissue_groups$kidney<-mp.df$kidsum
mp.tissue_groups$liver<-rowSums(mp.df[,grepl("liver", colnames(mp.df)), drop=FALSE])
mp.tissue_groups$pancreas<- rowSums(mp.df[,grepl("pancreas", colnames(mp.df)), drop=FALSE])
mp.tissue_groups$skin<-rowSums(mp.df[,grepl("skin", colnames(mp.df)), drop=FALSE])
mp.tissue_groups$spleen<-rowSums(mp.df[,grepl("spleen", colnames(mp.df)), drop=FALSE])
mp.tissue_groups$erythroblast<-rowSums(mp.df[,grepl("erythroblast", colnames(mp.df)), drop=FALSE])

        
m<-melt(aggregate(.~group, data=mp.tissue_groups[,-c(1)], mean), id.vars="group")
m$variable<-factor(m$variable, levels=c("macrophage", "dendritic", "eosonophil", "neutrophil", "monocyte", "NKCell", "TCell", "BCell", "erythroblast",
                                            "bladder", "kidney", "skin", "colon", "liver", "pancreas", "spleen"))

# FIGURES GROUP 1 ----
custom_palette=c("#9F66AB", #macrophage
                "#4582cc", #dendritic
                "black", #eosonophil
                "#8D782E", #neutrophil
                "#c96540", #monocyte
                "#D600C2", # NK
                "#3b312f", #TCell
                "#276647", #BCell
                 "red", #erythroblast
                "#6e655f", #bladder
                "#a11a2c", #kidney,
                "#66bf3d", #skin
                "#f57d7d", #colon
                "#4a1628", #liver
                "#00ccFF", #pancreas
                 "#2b3061") #spleen          

pdf(file="../Figures/tissue_comp.pdf",
    width=1.825, height=1.97, paper="special", bg="white", #used to be 1.75 x 1.78
    fonts="Helvetica", colormodel = "cmyk", pointsize=6)
ggplot(data=m)+
  geom_col(aes(x=group, y=value, fill=variable), position="stack", size=0.2, width=0.75)+
  n_samples(m, x=mp.tissue_groups$group, y=rep(0, nrow(mp.tissue_groups)), y_nudge=-0.03, fontsize=6)+
  geom_vline(xintercept = 5.5, linetype="dashed", size=0.25)+
  scale_fill_manual(values =custom_palette )+#ylim(c(0,1))+
  xlab(" ")+
  ylab("Tissue proportion")+
  theme_bw()+
  theme(plot.background=element_blank())+
  theme(#axis.title.x=element_blank(),
        legend.position = "none",
        axis.title=element_text(family="Helvetica", size=8),
        axis.text=element_text(family="Helvetica", size=6),
        axis.text.x = element_text(family="Helvetica", size=6, angle=25), #angle used to be 25 with hjust=1
        plot.title=element_text(family="Helvetica", size=6))
{dev.off()}

# FIGURES WITH TOTAL DNA ------------------------------------------------------------------------

#Of interest
biop<-droplevels(mp.df[grepl("BKV|Normal|Early", mp.df$group), ])
pdf(file="../Figures/Kidney_total_DNA_group_subset.pdf",
    width=1.5, height=1.5, paper="special", bg="white",
    fonts="Helvetica", colormodel = "cmyk", pointsize=6)
ggplot(data=biop, aes(x=group, y=kidsum*humanDNA*1000))+
  geom_boxplot(outlier.size = 0.5, color="black")+ geom_point(size=0.5)+
  n_samples(m, x=biop$group, y=rep(0, nrow(biop)), y_nudge=-2, fontsize=6)+
  geom_vline(xintercept = 3.5, linetype="dashed")+
  ylab("Kidney cfDNA (ng/ml)")+
  theme_bw()+coord_cartesian()+
  theme(plot.background=element_blank())+
  theme(legend.position="none")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(family="Helvetica", size=8),
        axis.text=element_text(family="Helvetica", size=6),
        axis.text.x=element_text(angle=25),
        plot.title=element_blank(), plot.margin = margin(t=5,r=1,l=5,b=5))
dev.off()

pdf(file="../Figures/WBC_total_DNA_group_subset.pdf",
    width=1.5, height=1.5, paper="special", bg="white",
    fonts="Helvetica", colormodel = "cmyk", pointsize=6)
ggplot(data=biop, aes(x=group, y=WBC*humanDNA*1000))+
  geom_boxplot(outlier.size = 0.5, color="black")+ geom_point(size=0.5)+
  n_samples(m, x=biop$group, y=rep(0, nrow(biop)), y_nudge=-9, fontsize=6)+
  geom_vline(xintercept = 3.5, linetype="dashed")+
  ylab("Leukocyte cfDNA (ng/ml)")+
  theme_bw()+coord_cartesian()+
  theme(plot.background=element_blank())+
  theme(legend.position="none")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(family="Helvetica", size=8),
        axis.text=element_text(family="Helvetica", size=6),
        axis.text.x=element_text(angle=25),
        plot.title=element_blank(), plot.margin = margin(t=5,r=1,l=5,b=5))
dev.off()

# STATISTICS ----
BKVN<-mp.df[mp.df$group=="BKVN", ]
BKVNneg<-mp.df[mp.df$group=="BKV+/N-", ]
NORMAL<-mp.df[mp.df$group=="Normal", ]
NOUTI<-mp.df[mp.df$group=="no-UTI", ]
UTI<-mp.df[mp.df$group=="UTI", ]
EARLY<-mp.df[mp.df$group=="Early", ]

mean(NORMAL$WBC*NORMAL$humanDNA)*1000
#BKVN vs Early
z<-rbind(BKVN, EARLY)
wilcox.test(z$kidsum*z$humanDNA~z$group)
wilcox.test(z$WBC*z$humanDNA~z$group)
wilcox.test(z$bladder*z$humanDNA~z$group)

#BKV+/N- vs EARLY
y<-rbind(BKVNneg, EARLY)
wilcox.test(y$kidsum*y$humanDNA~y$group)
wilcox.test(y$WBC*y$humanDNA~y$group)
wilcox.test(y$bladder*y$humanDNA~y$group)

# BKVN vs BKV+/N-
a<-rbind(BKVN, BKVNneg)
wilcox.test(a$kidsum*a$humanDNA~a$group)
wilcox.test(a$WBC*a$humanDNA~a$group)
wilcox.test(a$bladder*a$humanDNA~a$group)

# BKVN vs Normal
b<-rbind(BKVN, NORMAL)
wilcox.test(b$kidsum*b$humanDNA~b$group)
wilcox.test(b$WBC*b$humanDNA~b$group)
wilcox.test(b$bladder*b$humanDNA~b$group)

# BKV+/N- vs Normal
c<-rbind(BKVNneg, NORMAL)
wilcox.test(c$kidsum*c$humanDNA~c$group)
wilcox.test(c$WBC*c$humanDNA~c$group)
wilcox.test(c$bladder*c$humanDNA~c$group)

#UTI vs No UTI
d<-rbind(UTI, NOUTI)
wilcox.test(d$kidsum*d$humanDNA~d$group)
wilcox.test(d$WBC*d$humanDNA~d$group)
wilcox.test(d$bladder*d$humanDNA~d$group)

#UTI vs Normal
e<-rbind(UTI, NORMAL)
wilcox.test(e$kidsum*e$humanDNA~e$group)
wilcox.test(e$WBC*e$humanDNA~e$group)
wilcox.test(e$bladder*e$humanDNA~e$group)

mean(UTI$WBC*UTI$humanDNA*1000)
mean(NORMAL$WBC*NORMAL$humanDNA*1000)
#No UTI vs NORMAL
f<-rbind(NOUTI, NORMAL)
wilcox.test(f$kidsum*f$humanDNA~f$group)
wilcox.test(f$WBC*f$humanDNA~f$group)
wilcox.test(f$bladder*f$humanDNA~f$group)

# BKVN vs UTI
g<-rbind(BKVN, UTI)
wilcox.test(g$kidsum~g$group)
