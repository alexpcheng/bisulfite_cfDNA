library(data.table)
mp.df<-fread("../V2/tissues_of_origin/mp.txt")
mp2<-data.frame(t(mp.df[, 2:114]))
mp2$thing<-rownames(mp2)
mp2$thing<-gsub('[0-9]+', '', mp2$thing)
mp3<-aggregate(.~thing, mp2, sum)
mp4<-data.frame(t(mp3))
colnames(mp4)<-mp3$thing
mp5<-mp4[-1,]
mp5$meth_id<-mp.df$sample

sample_groups<-fread("../lists/sample_groups.txt", header=FALSE)
colnames(sample_groups)<-c("meth_id", "group")
mp6<-merge(mp5, sample_groups, by="meth_id")
mp6$group<-gsub("HLY", "no UTI", mp6$group)
mp6$group<-gsub("vir-/inf-", "normal biopsy", mp6$group)

sample_names<-fread("../../KTxBS_renamed_for_safety/lists/sample_names.txt")
colnames(sample_names)[1]<-"meth_id"
mp7<-merge(sample_names, mp6, by="meth_id")

microbial_fraction<-fread("../V1/read_statistics/BS_treated_samples.txt")
microbial_fraction$human_fraction<-1-(microbial_fraction$microbe_mapped/microbial_fraction$total_reads)
microbial_fraction<-microbial_fraction[, c("sample", "human_fraction")]
colnames(microbial_fraction)[1]<-"meth_id"

cfDNA<-fread("../lists/sample_cfDNA_extracts")
cfDNA<-merge(cfDNA, microbe_fraction, by="meth_id")
cfDNA$humanDNA<-cfDNA$cfDNA_ng_ul*cfDNA$elution_volume/cfDNA$urine_volume_ul #*cfDNA$human_fraction

mp8<-merge(mp7, cfDNA[, c("meth_id", "humanDNA")], by="meth_id")

fwrite(file = "../V1/mp_agg.txt", x = mp8, quote = FALSE, sep = '\t', row.names = FALSE)