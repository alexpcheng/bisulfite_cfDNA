### Burnham, De Vlaminck: 2017
# donorfraction.autosomalandY.R

#Slight modifications from original script done by Alexandre Pellan Cheng
#Adapted for BS treated samples
### libraries
suppressMessages(library(HMMcopy))

### Variables
args <- commandArgs(trailingOnly = T) # command line arguments
sam = as.character(args[1])
sample = paste0("../V1/donor_fractions/",sam)

#import list of removal bins for Y-chromosome
removals <- as.numeric(as.matrix(read.table("../lists/BKVcohort_ff_removal.list",header = F, sep="\t")))

########################################################
#known F->F from UTI cohort

# Function -----------------------------------------------------------------------------------------------------------

get_sex_frac <- function(sample,thr=0.8, stdv_multi = 6){
  
  #import files, references, and execute analysis using HMMcopy format and procedure
  rfile <- paste0(sample,".auto.readcounts.wig")
  gfile <- "../Bin/HMMcopy/data/500.hg19.meth.gc.wig"
  mfile <- "../Bin/HMMcopy/data/500.hg19.full.map.wig"
  normal_reads <- wigsToRangedData(readfile = rfile, gcfile = gfile,mapfile = mfile, verbose=F)
  normal_copy <- correctReadcount(normal_reads,verbose = F, samplesize = 50000)
  
  #change to easy to read data.frame format
  df <- cbind(data.frame(normal_copy@ranges),normal_copy$cor.map, normal_copy$map)
  df <- df[,c(2,3,6,7)]
  colnames(df)<-c("chrom","binstart","value","mapa")
  df[is.na(df)] <- 0
  
  #determine cutoffs to eliminate overrepresented bins
  ythresh <- mean(df[df$chrom=="chrY"&!is.na(df$value),]$value) + stdv_multi*sd(df[df$chrom=="chrY"&!is.na(df$value),]$value)
  athresh <- mean(df[df$chrom!="chrX"&df$chrom!="chrY"&!is.na(df$value),]$value) + stdv_multi*sd(df[df$chrom!="chrX"&df$chrom!="chrY"&!is.na(df$value),]$value)
  
  #determine the sub.frames for each sex chromosome as well as automsomes
  ysa <- df[df$chrom=="chrY" & df$mapa > thr & df$value <ythresh &!( df$binstart %in% removals),]
  asa <- df[df$chrom!="chrX" & df$chrom!="chrY" & df$mapa > thr & df$value <athresh,]
 
  #calculate donor fraction based on the relative coverages of chrY and the average autosomal coverage.
  ratta.M <- (2*100*mean(ysa$value)/mean(asa$value))
  ratta.F <- 100*(1-(2*mean(ysa$value)/mean(asa$value)))
  return(c(ratta.M, ratta.F))
}

## Code execution #################################################

result = get_sex_frac(sample=sample)
names <- c("Donor_Gender","Donor_Fraction")
report.M <- cbind(names, c("M",as.numeric(result[1])))
report.F <- cbind(names, c("F",as.numeric(result[2])))

#export file
dir.create("../V2/donor_fractions", showWarnings=F)
sample_out<-paste0("../V2/donor_fractions/", sam)
write.table(report.M, paste0(sample_out,".M.autoY.sexmm"),quote = F,row.names = F, col.names = F, sep="\t")
write.table(report.F, paste0(sample_out,".F.autoY.sexmm"),quote = F,row.names = F, col.names = F, sep="\t")
