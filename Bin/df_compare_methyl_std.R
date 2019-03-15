#!/usr/bin/env RScript
# Title: Method comparison for donor fractions
# Authors: Alexandre Pellan Cheng
# Brief description: Compares donor fraction measurements obtained via standard sequencing and WGBS

# Creates supplementary figure [X]
# Requires succesful running of donor_fractions_make_list.R (you must first have calculated all donor fractions)
# Assumes files is run from project_name/Bin

# Initialize -------------------------------------------------------------------------------------------------
rm(list=ls())

# Load libraries ---------------------------------------------------------------------------------------------
library(data.table)
library(ggplot2)

# Paths ------------------------------------------------------------------------------------------------------
lists <- "../lists/"
donor_fractions_path<-"../V2/donor_fractions/"
meth_donor_fraction <- "../V2/donor_fractions/donor_fractions.txt"

# Load data --------------------------------------------------------------------------------------------------
donor_fractions_meth<-fread(meth_donor_fraction)
# See Burnham 2018 (Nature Communications) for methods on obtaining donor fractions from non sex mismatched
donor_fractions_std<-fread(paste0(lists, "donor_fractions_std.list"))
colnames(donor_fractions_std)<-c("sample", "std_df")
donor_fractions<-merge(donor_fractions_meth, donor_fractions_std, by="sample")

save_eps=TRUE
core<-cor.test(donor_fractions$std_df, as.numeric(donor_fractions$donor_fractions), method="spearman")
if(save_eps){pdf(file="../Figures/DF_comparison.pdf",
                 width=3, height=3, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=8)}
ggplot(donor_fractions) + 
  geom_point(aes(x=std_df/100, y=as.numeric(donor_fractions)/100))+
  geom_segment(aes(x=0, y=0, xend=100, yend=100), linetype="dashed")+
  xlab("Standard sequencing donor fraction")+
  ylab("WGBS donor fraction (%)")+
  theme_bw()+coord_cartesian(ylim=c(0,1), xlim=c(0,1))+
  theme(plot.background=element_blank())+
  theme(legend.position="none")+
  theme(axis.title.x=element_text(family="Helvetica", size=8), axis.title.y=element_text(family="Helvetica", size=8),
        axis.text=element_text(family="Helvetica", size=8), axis.text.x=element_text(angle=0),
        plot.title=element_text(family="Helvetica", size=8))
if(save_eps){dev.off()}

