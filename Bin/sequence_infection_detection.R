# Load libraries ----------------------------------------------------------------------------------------------------------
library(data.table)
library(plyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
# Paths -------------------------------------------------------------------------------------------------------------------
lists<- "../lists/"
infection.path<-"../V1/infections/"
mixing_parameters.path<-"../V2/tissues_of_origin/"

# Functions ---------------------------------------------------------------------------------------------------------------
get_infection<-function(filename, infection.path, common_microbe_contaminants){
  file<-data.frame(fread(paste0(infection.path, filename, ".grammy.tab"), header=TRUE))
  f=filename
  species.abundance<-file[,c("superkingdom", "species", "RelCoverage", "genus")]
  
  species.abundance<-merge(species.abundance, common_microbe_contaminants, by="genus", all.x=TRUE)
  species.abundance<-species.abundance[is.na(species.abundance$genus_name), ]
  species.abundance<-species.abundance[, c("superkingdom", "species", "RelCoverage")]
  
  colnames(species.abundance)[3]<-"RelCoverage_meth_seq"
  species.abundance<-aggregate(.~species+superkingdom, data=species.abundance, sum)
  species.abundance$meth_id<-f
  species.abundance<-species.abundance[(species.abundance$superkingdom!=2759), ]
  species.abundance<-species.abundance[(species.abundance$superkingdom!=2157), ]
  
  bact<-species.abundance[species.abundance$superkingdom==2, ]
  viru<-species.abundance[species.abundance$superkingdom==10239, ]
  species.abundance<-rbind(bact[which.max(bact$RelCoverage_meth_seq), ],
                           viru[which.max(viru$RelCoverage_meth_seq), ])
  
  return(species.abundance)
}

# Load samples -------------------------------------------------------------------------------------------------------------------
sams<-fread(paste0(lists, "included_samples.txt"), header=F, col.names = "meth_id")
sample_groups<-fread("../lists/sample_groups.txt", header=FALSE, col.names = c("meth_id", "group"))
common_microbe_contaminants<-fread(paste0(lists, "common_sequence_contaminants.toremove.txt"), header=TRUE, col.names=c("genus", "genus_name"))
mp.df<-fread(paste0(mixing_parameters.path, "mp.txt"))
colnames(mp.df)[1]<-"meth_id"

wgbs.rge<-rbindlist(lapply(X = as.list(sams$meth_id), FUN = get_infection, infection.path, common_microbe_contaminants))

#Reformat wgbs.rge
wgbs.rge<-merge(wgbs.rge[wgbs.rge$superkingdom==2, ], wgbs.rge[wgbs.rge$superkingdom==10239, ], by="meth_id")
colnames(wgbs.rge)<-c("meth_id","bacteria_species", "bacteria_superkingdom", "bacteria_RGE", "virus_species", "virus_superkingdom", "virus_RGE")

#wgbs.rge<-wgbs.rge[!grepl("MET-1-52|MET-1-56|MET-1-57|MET-1-59|MET-1-60|MET-1-63", wgbs.rge$meth_id), ]

wgbs.rge$potential_uti<-"-NO UTI-"
wgbs.rge$potential_uti[wgbs.rge$bacteria_RGE>0.08975]<-"-UTI-"
wgbs.rge$potential_virus<-"-NO VIRUS-"
wgbs.rge$potential_virus[wgbs.rge$virus_RGE>10**3]<-"-VIRUS-"
tt<-merge(wgbs.rge, sample_groups, by="meth_id")
wgbs.rge$label<-paste0(wgbs.rge$potential_uti, wgbs.rge$potential_virus, wgbs.rge$meth_id)

infection_samples<-merge(mp.df, wgbs.rge[, c("meth_id", "label")], by="meth_id")
infection_samples<-infection_samples[!(grepl("-NO UTI--NO VIRUS", infection_samples$label)), ]
infection_samples<-infection_samples[, c(1:91, 93:99, 101:115, 123)] #only grab primary references and label column

col.order<-hclust(dist(infection_samples[,c(2:113)], method = "man"), method="complete")
infection_samples<-infection_samples[col.order$order, ]

m<-melt(infection_samples[, -c("meth_id")], id.vars="label")

summed_tissues<-fread("../V1/mp_agg.txt")
summed_tissues<-merge(summed_tissues, infection_samples[,c("meth_id")], by="meth_id")
summed_tissues<-summed_tissues[, c(6:21)]
colnames(summed_tissues)[order(colMeans(summed_tissues))]
library(ggdendro)

dendr<-dendro_data(col.order, type="rectangle")

pdf(file="../Figures/infection_determined_seq_dendro.pdf",
    width=0.75, height=2.7, paper="special", bg="white",
    fonts="Helvetica", colormodel = "cmyk", pointsize=8)
ggplot()+geom_segment(data=segment(dendr), aes(x=y, y=x, xend=yend, yend=xend))+
  theme(plot.background=element_blank(),
          legend.position="none",
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank())
dev.off()

m$label<-factor(m$label, levels=infection_samples$label)
m$color_label<-gsub("MET.*", "", m$label)
m$placeholder="placeholder"

m$cell_type<-gsub("[0-9]+", "", m$variable)
m2<-m[, c("label", "value", "color_label", "placeholder", "cell_type")]
m2<-aggregate(.~(label+cell_type+placeholder+color_label), m2, sum)
#m2$cell_type<-factor(m2$cell_type, row.order$labels)

#Figure ----
pdf(file="../Figures/infection_determined_seq_too1.pdf",
    width=1.5, height=2.7, paper="special", bg="white",
    fonts="Helvetica", colormodel = "cmyk", pointsize=8)
std_theme<-theme(plot.background=element_blank(),
                 legend.position="none",
                 #axis.text = element_text(angle=90),
                 axis.text = element_blank(),
                 axis.title = element_blank(),
                 axis.ticks = element_blank(),
                 panel.grid = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_blank())
ggplot(data=m2)+
  geom_tile(aes(y=label, x=cell_type, fill=value), color="black")+
  scale_fill_distiller(palette="RdYlBu")+
  #scale_fill_gradient2(low="purple", mid="white", high="green", midpoint=0.5, breaks=c(0,0.5,1))+
  std_theme
dev.off()
max(m2$value)
min(m2$value)
pdf(file="../Figures/infection_determined_seq_too2.pdf",
    width=0.3, height=2.7, paper="special", bg="white",
    fonts="Helvetica", colormodel = "cmyk", pointsize=8)
std_theme<-theme(plot.background=element_blank(),
                 legend.position="none",
                 axis.text = element_blank(),
                 axis.title = element_blank(),
                 axis.ticks = element_blank(),
                 panel.grid = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_blank())
ggplot(data=m2)+geom_tile(aes(y=label, x=placeholder, fill=color_label), color="black", size=0.2)+
  scale_fill_manual(values=c("red", "#0080ff", "green"))+
  std_theme
dev.off()

