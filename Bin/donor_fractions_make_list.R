#!/usr/bin/env RScript
# Title: Donor fraction Get
# Authors: Alexandre Pellan Cheng
# Brief description: Get the donor fraction from each sex mismatched patient
# Called by donor_fractions_par.sh

# Load libraries ---------------------------------------------------------------------------------------------
library(data.table)
# Paths ------------------------------------------------------------------------------------------------------
donor_genders.list <- "../lists/only_sexmm.list"
donor_fractions.path <- "../V2/donor_fractions/"
# Load data --------------------------------------------------------------------------------------------------
donor_genders<-fread(donor_genders.list, header=F)
colnames(donor_genders)<-c("sample", "donor_gender", "recipient_gender")

# Create dataframe with sample and donor fraction ------------------------------------------------------------
donor_fractions<-data.frame(matrix(ncol=2, nrow=nrow(donor_genders)))
colnames(donor_fractions)<-c("sample", "donor_fractions")

for (i in 1:nrow(donor_genders)){
  file<-paste0("../V2/donor_fractions/",
               donor_genders$sample[i], ".",
               donor_genders$donor_gender[i],
               ".autoY.sexmm")
  df<-fread(file, header = F, colClasses = c("character"))
  sample<-donor_genders$sample[i]
  fraction<-df$V2[2]
  donor_fractions[i,1]<-sample
  donor_fractions[i,2]<-fraction
}

# Write data -------------------------------------------------------------------------------------------------
fwrite(donor_fractions,
       file=paste0(donor_fractions.path, "donor_fractions.txt"),
       sep="\t")
