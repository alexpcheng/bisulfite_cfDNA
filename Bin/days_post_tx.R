#!/usr/bin/env RScript
# Title: Tissue composition by pathology
# Authors: Alexandre Pellan Cheng
# Brief description: Creates figure XYZ

# Initialize -------------------------------------------------------------------------------------------------
rm(list=ls())

# Load libraries ---------------------------------------------------------------------------------------------
library(ggplot2)
library(data.table)

# Load files -------------------------------------------------------------------------------------------------

days_since_tx<-fread("../lists/days_since_tx", header=TRUE)
days_since_tx$DAYS_POST_TRANSPLANT<-as.numeric(days_since_tx$DAYS_POST_TRANSPLANT)

sample_groups<-fread("../lists/sample_groups.txt", header=FALSE, col.names=c("sample", "group"))

mp.df<-fread("../V2/tissues_of_origin/mp.txt")

cfDNA<-fread("../lists/sample_cfDNA_extracts")
colnames(cfDNA)[1]<-"sample"
microbial_fraction<-fread("../V1/read_statistics/BS_treated_samples.txt")
microbial_fraction$human_fraction<-1-(microbial_fraction$microbe_mapped/microbial_fraction$total_reads)
microbial_fraction<-microbial_fraction[, c("sample", "human_fraction")]

df<-Reduce(f = merge, x = list(days_since_tx, sample_groups, mp.df, cfDNA, microbial_fraction))

df$humanDNA<-df$cfDNA_ng_ul*df$elution_volume/df$urine_volume_ul*df$human_fraction

df<-df[, c("sample", "kidsum", "humanDNA", "group.x", "DAYS_POST_TRANSPLANT")]

df<-df[grepl("ETP|vir-|HLY", df$group.x), ]

summed<-fread("../../mp_agg.txt")
colnames(summed)[1]<-"sample"
df<-merge(df[, c("DAYS_POST_TRANSPLANT", "sample")], summed, by="sample")

df$macrophage<-df$macrophage*df$humanDNA
df$dendritic<-df$dendritic*df$humanDNA
df$eosonophil<-df$eosonophil*df$humanDNA
df$neutrophil<-df$neutrophil*df$humanDNA
df$monocyte<-df$monocyte*df$humanDNA
df$NK<-df$NK*df$humanDNA
df$TCell<-df$TCell*df$humanDNA
df$BCell<-df$BCell*df$humanDNA
df$erythroblast<-df$erythroblast*df$humanDNA
df$bladder<-df$bladder*df$humanDNA
df$kidney<-df$kidney*df$humanDNA
df$skin<-df$skin*df$humanDNA
df$colon<-df$colon*df$humanDNA
df$liver<-df$liver*df$humanDNA
df$pancreas<-df$pancreas*df$humanDNA
df$spleen<-df$spleen*df$humanDNA

df2<-df[, c(2,7:22)]
df2<-aggregate(.~DAYS_POST_TRANSPLANT, df2, mean)

occurences<-df[ ,c("DAYS_POST_TRANSPLANT"), drop=FALSE]
occurences$n<-1
occurences<-aggregate(.~DAYS_POST_TRANSPLANT, occurences, sum)

occurences$n[occurences$n>1]<-2

m<-melt(df2, id.vars="DAYS_POST_TRANSPLANT")
m2<-aggregate(.~DAYS_POST_TRANSPLANT, m[, c("DAYS_POST_TRANSPLANT", "value")], sum)
m2<-merge(m2, occurences, by="DAYS_POST_TRANSPLANT")
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
m$variable<-factor(m$variable, levels=c("macrophage", "dendritic", "eosonophil", "neutrophil", "monocyte", "NKCell", "TCell", "BCell", "erythroblast",
                                        "bladder", "kidney", "skin", "colon", "liver", "pancreas", "spleen"))

pdf(file="../Figures/tissue_comp_by_days.pdf",
    width=2, height=1.5, paper="special", bg="white",
    fonts="Helvetica", colormodel = "cmyk", pointsize=6)
ggplot()+
  geom_area(data=m, aes(x=DAYS_POST_TRANSPLANT, y=value*1000, fill=variable), position='stack')+
  geom_point(data=m2, aes(x=DAYS_POST_TRANSPLANT, y=value*1000, shape =as.factor(n)), size=1)+
  geom_vline(xintercept=5, linetype="dashed", size=0.5)+
  geom_vline(xintercept=100, linetype="dashed", size=0.5)+
  scale_x_log10()+
  scale_fill_manual(values=custom_palette)+
  xlab("Days post transplant")+
  ylab("Tissue cfDNA (ng/ml)")+#ylim(c(0,1))+
  theme_bw()+
  theme(plot.background=element_blank())+
  theme(legend.position = "none",
        axis.title=element_text(family="Helvetica", size=8),
        #axis.title.y=element_blank(),
        axis.text=element_text(family="Helvetica", size=6),
        axis.text.x = element_text(family="Helvetica", size=6),
        plot.title=element_blank(), plot.margin = margin(t=5,r=1,l=5,b=5))
dev.off()
 
