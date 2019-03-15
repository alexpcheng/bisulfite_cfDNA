#!/usr/bin/env RScript
# Title: Unsupervised clustering of reference methylomes and sample projection
# Authors: Alexandre Pellan Cheng
# Brief description: Generates figures X and Y for [publication]

# Assume setwd = project_folder/Bin
#Initialize console ----------------------------------------------------------
rm(list=ls())

#libraries --------------------------------------------------------------------
library(data.table)
library(parallel)
library(ggplot2)
library(umap)
library(dplyr)
library(ggpubr)
#Paths ------------------------------------------------------------------------

#references.path <-"../Methylation_References/"
references.path <- "../Methylation_References_snaked/MethylMatrix/"
samples.path<-"../V1/binned_samples/"

#Functions --------------------------------------------------------------------
# Creates a list with each element containing the reference matrix or a sample.
  # NANs are removed from dataframes and sample methylation percentages are 
  # scaled to a fraction to be consistent with the reference matrix.
process_files<-function(filename, samples.path, references.path){
  if(filename=="MethylMatrix_binned"){
    col<-fread("../Methylation_References_snaked/lookup_table.txt", header=F)
    tiled<-fread(paste0(references.path, filename))
    tiled<-tiled[complete.cases(tiled),]
    colnames(tiled)<-c("chr", "start", "end", col$V4)
  }
  else if(grepl("MET", filename)){
    tiled<-fread(paste0(samples.path, filename))
    colnames(tiled)<-c("chr", "start", "end", "pmeth")
    tiled<-tiled[complete.cases(tiled),] #Remove NAN
    tiled$pmeth<-tiled$pmeth/100 #Normalize to 1
    colnames(tiled)[4]<-filename
  }
  return(tiled)
}

# Renames tissues for easier plotting
tissue_names<-function(sample_row, coll){
  if (grepl("kidney", sample_row[coll]))
    colr<-"kidney"
  else if (grepl("mesangial", sample_row[coll]))
    colr<-"kidney_mesangial"
  else if (grepl("podocyte", sample_row[coll]))
    colr<-"kidney_podocyte"
  else if (grepl("Large_intestine|sigmoid", sample_row[coll]))
    colr<-"gut_large_intestine"
  else if (grepl("small_intestine", sample_row[coll]))
    colr<-"gut_small_intestine"
  else if (grepl("bladder", sample_row[coll]))
    colr<-"zbladder" #for layering
  else if (grepl("skin", sample_row[coll]))
    colr<-"skin"
  else if (grepl("liver", sample_row[coll]))
    colr<-"liver"
  else if (grepl("hepatocyte", sample_row[coll]))
    colr<-"liver_hepatocyte"
  else if (grepl("pancreas", sample_row[coll]))
    colr<-"pancreas_z"
  else if (grepl("islet", sample_row[coll]))
    colr<-"pancreas_islet"
  else if (grepl("MET", sample_row[coll]))
    colr<-"urine sample"
  else if (grepl("macrophage", sample_row[coll]))
    colr<-"macrophage"
  else if (grepl("monocyte", sample_row[coll]))
    colr<-"monocyte"
  else if (grepl("dendritic", sample_row[coll]))
    colr<-"dendritic"
  else if (grepl("eosonophil", sample_row[coll]))
    colr<-"eosonophil"
  else if (grepl("neutrophil", sample_row[coll]))
    colr<-"neutrophil"
  else if (grepl("NK", sample_row[coll]))
    colr<-"zNKell"
  else if (grepl("TCell", sample_row[coll]))
    colr<-"Tcell"
  else if (grepl("BCell", sample_row[coll]))
    colr<-"Bcell"
  else if (grepl("spleen", sample_row[coll]))
    colr<-"spleen"
  else if (grepl("erythroblast", sample_row[coll]))
    colr<-"erythroblast"
  else
    colr<-sample_row[coll]
  return(colr)
}

#Create sample list
sample_list<-list.files(path = samples.path, pattern="*MET")
#references_list<-list.files(path=references.path, pattern="MM.binned")
references_list<-list.files(path=references.path, pattern="MethylMatrix_binned")
samples_and_references.locations<-c(references_list, sample_list)

#Select common regions to references and samples ------------------------------
refs_and_sams.list<-mclapply(samples_and_references.locations, FUN=process_files, samples.path, references.path, mc.cores=10)
refs_and_sams.df<-Reduce(function(x, y) merge(x,y, by=c("chr", "start", "end")), refs_and_sams.list)

refs_and_sams.features<-refs_and_sams.df[,4:(ncol(refs_and_sams.df))] #removes chr/start/end
num_samples<-length(sample_list)
num_references<-ncol(refs_and_sams.features)-num_samples
refs.features<-refs_and_sams.features[,1:num_references]

# Perform unsupervised clustering on references -------------------------------
# Kmeans
set.seed(1)
km<-kmeans(t(refs.features), centers=4, nstart = 10)

#PCA
pca<-prcomp(t(refs.features))
summary(pca) # Get % variance for PC1 and P
pca_components<-as.data.frame(pca$x)
pca_components$sample<-factor(rownames(pca_components), levels=rownames(pca_components))
pca_components$tissue<-apply(X = pca_components, MARGIN = 1, FUN = tissue_names, ncol(pca_components))
pca_components$cluster<-km$cluster
pca_components<-pca_components[order(km$cluster),]

# UMAP
UMAP=umap(t(refs.features), random_state=3) #random state set for consistency
UMAP_dims<-data.frame(UMAP$layout)
UMAP_dims$sample<-factor(colnames(refs.features), levels=colnames(refs.features))
UMAP_dims$tissue<-apply(X = UMAP_dims, MARGIN = 1, FUN = tissue_names, ncol(UMAP_dims))
UMAP_dims$cluster<-km$cluster

# PLOT REFERENCES
custom_palette=c("#276647", #BCell
                 "#4582cc", #dendritic
                 "black", #eosonophil
                 "red", #erythroblast
                 "#f57d7d", "#990000", #large and small intestine
                 "#a11a2c", "#a14c1a", "#003333", #kidney, mesangial, podocyte
                 "#4a1628", "#0000cc", #liver, hepatocyte
                 "#9F66AB", #"#99CC00", #macrophage
                 "#c96540", #monocyte
                 "#8D782E", #"#d4be42", #neutrophil
                 "#0EA08E", "#00ccFF", #islet (#009999), pancreas
                 "#66bf3d", #skin
                 "#2b3061", #spleen
                 "#3b312f", #TCell
                 "#6e655f", #bladder
                 "#D600C2")#  "#f0c295") #NK                 
cluster_palette=c("#574AE2", "#00a5CF", "#DE1A1A", "#29BF12")

pdf(file="../Figures/PCA_refs.pdf",
    width=3, height=3, paper="special", bg="white",
    fonts="Helvetica", colormodel = "cmyk", pointsize=6)
ggplot(data=pca_components)+
  stat_ellipse(geom="polygon", alpha=0.5, aes(x=PC1, y=PC2, fill=factor(cluster)))+scale_fill_manual(values=cluster_palette)+
  geom_point(aes(x=PC1, y=PC2, color=tissue))+scale_color_manual(values=custom_palette)+
  theme_bw()+
  xlab("PC1 (29.9%)")+
  ylab("PC2 (16.8%)")+
  theme(plot.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.title=element_text(family="Helvetica", size=8), 
        axis.text=element_text(family="Helvetica", size=6),
        plot.title=element_text(family="Helvetica", size=6))+
  theme(legend.position = "none")
dev.off()

pdf(file="../Figures/UMAP_refs.pdf",
    width=3, height=3, paper="special", bg="white",
    fonts="Helvetica", colormodel = "cmyk", pointsize=6)
ggplot(data=UMAP_dims %>% arrange(tissue))+
  stat_ellipse(geom="polygon",alpha=0.5, aes(x=X1, y=X2, fill=factor(cluster)))+scale_fill_manual(values=cluster_palette)+
  geom_point(aes(x=X1, y=X2, color=tissue), size=1)+scale_color_manual(values=custom_palette)+
  theme_bw()+xlab("UMAP1")+ylab("UMAP2")+
  theme(plot.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.title=element_text(family="Helvetica", size=8),
        axis.text=element_text(family="Helvetica", size=6),
        plot.title=element_text(family="Helvetica", size=6))+
 theme(legend.position = "none")
dev.off()

#Project data onto maps
df<-data.frame(fread("../V2/donor_fractions/donor_fractions.txt"))
df<-df[!(grepl("MET-1-10|MET-1-18", df$sample)),] #remove samples that also have a bone marrow transplant

sams.features<-refs_and_sams.features[,(num_references+1):(ncol(refs_and_sams.features))]

#PCA
#sam_on_PCA<-data.frame(scale(t(sams.features), pca$center, pca$scale)%*% pca$rotation)
#sam_on_PCA$sample<-rownames(sam_on_PCA)
#sam_on_PCA<-merge(sam_on_PCA, df, by="sample")

sam_on_PCA<- data.frame(predict(pca, t(sams.features)))
sam_on_PCA$sample<-colnames(sams.features)
sam_on_PCA<-merge(sam_on_PCA, df, by="sample")


#UMAP
sample.umap = predict(UMAP, t(sams.features))

sam_on_UMAP<-data.frame(sample.umap)
sam_on_UMAP$sample<-colnames(sams.features)
sam_on_UMAP<-merge(sam_on_UMAP, df, by="sample")

pdf(file="../Figures/PCA_refs_sams_DF.pdf",
    width=1.825, height=1.825, paper="special", bg="white", #used to be 1.75 by 1.75
    fonts="Helvetica", colormodel = "cmyk")
ggplot(data=pca_components)+
  stat_ellipse(geom="polygon", alpha=0.5, aes(x=PC1, y=PC2, fill=factor(cluster)))+scale_fill_manual(values=cluster_palette)+
  theme_bw()+
  xlab("PC1 (29.9%)")+
  ylab("PC2 (16.8%)")+
  theme(plot.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(family="Helvetica", size=8, margin=margin(t=0)), 
        axis.title.y=element_text(family="Helvetica", size=8, margin=margin(r=0)),
        axis.text=element_text(family="Helvetica", size=6),
        plot.title=element_text(family="Helvetica", size=6), legend.position = "none")+
  geom_point(data=sam_on_PCA, aes(x=PC1, y=PC2, color=donor_fractions), size=1)+scale_color_gradient(low="white", high="black")+
  geom_point(data=sam_on_PCA, aes(x=PC1, y=PC2), fill="NA", color="black", pch=21, size=1.5, stroke=0.1)
dev.off()

pdf(file="../Figures/UMAP_refs_sams_DF.pdf",
    width=1.825, height=1.825, paper="special", bg="white", #used to be 1.75 by 1.75
    fonts="Helvetica", colormodel = "cmyk")
ggplot(data=UMAP_dims)+
  stat_ellipse(geom="polygon", alpha=0.5, aes(x=X1, y=X2, fill=factor(cluster)))+scale_fill_manual(values=cluster_palette)+
  theme_bw()+xlab("UMAP1")+ylab("UMAP2")+
  theme(plot.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(family="Helvetica", size=8, margin=margin(t=0)), 
        axis.title.y=element_text(family="Helvetica", size=8, margin=margin(r=0)),
        axis.text=element_text(family="Helvetica", size=6),
        plot.title=element_text(family="Helvetica", size=6),
        legend.position = "none")+
  geom_point(data=sam_on_UMAP, aes(x=X1, y=X2, color=donor_fractions), size=1)+scale_color_gradient(low="white", high="black")+
  geom_point(data=sam_on_UMAP, aes(x=X1, y=X2), fill="NA", color="black", pch=21, size=1.5, stroke=0.1)
dev.off()

# Additional principal components ---------------------------------------------------
std_theme<-theme(plot.background=element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.title.x=element_text(family="Helvetica", size=8), 
                axis.title.y=element_text(family="Helvetica", size=8),
                axis.text=element_text(family="Helvetica", size=6),
                plot.title=element_text(family="Helvetica", size=6),
                legend.position = "none")

g12<-ggplot(data=pca_components)+
  stat_ellipse(geom="polygon", alpha=0.5, aes(x=PC1, y=PC2, fill=factor(cluster)))+scale_fill_manual(values=cluster_palette)+
  geom_point(aes(x=PC1, y=PC2, color=tissue))+scale_color_manual(values=custom_palette)+
  theme_bw()+
  xlab("PC1 (29.9%)")+
  ylab("PC2 (16.8%)")+
  std_theme

g21<-ggplot(data=pca_components)+
  stat_ellipse(geom="polygon", alpha=0.5, aes(y=PC1, x=PC2, fill=factor(cluster)))+scale_fill_manual(values=cluster_palette)+
  geom_point(aes(y=PC1, x=PC2, color=tissue))+scale_color_manual(values=custom_palette)+
  theme_bw()+
  ylab("PC1 (29.9%)")+
  xlab("PC2 (16.8%)")+
  std_theme

g31<-ggplot(data=pca_components)+
  stat_ellipse(geom="polygon", alpha=0.5, aes(y=PC1, x=PC3, fill=factor(cluster)))+scale_fill_manual(values=cluster_palette)+
  geom_point(aes(y=PC1, x=PC3, color=tissue))+scale_color_manual(values=custom_palette)+
  theme_bw()+
  ylab("PC1 (29.9%)")+
  xlab("PC3 (11.9%)")+
  std_theme

g13<-ggplot(data=pca_components)+
  stat_ellipse(geom="polygon", alpha=0.5, aes(x=PC1, y=PC3, fill=factor(cluster)))+scale_fill_manual(values=cluster_palette)+
  geom_point(aes(x=PC1, y=PC3, color=tissue))+scale_color_manual(values=custom_palette)+
  theme_bw()+
  xlab("PC1 (29.9%)")+
  ylab("PC3 (11.9%)")+
  std_theme

g23<-ggplot(data=pca_components)+
  stat_ellipse(geom="polygon", alpha=0.5, aes(x=PC2, y=PC3, fill=factor(cluster)))+scale_fill_manual(values=cluster_palette)+
  geom_point(aes(x=PC2, y=PC3, color=tissue))+scale_color_manual(values=custom_palette)+
  theme_bw()+
  xlab("PC2 (16.8%)")+
  ylab("PC3 (11.9%)")+
  std_theme

g32<-ggplot(data=pca_components)+
  stat_ellipse(geom="polygon", alpha=0.5, aes(x=PC3, y=PC2, fill=factor(cluster)))+scale_fill_manual(values=cluster_palette)+
  geom_point(aes(x=PC3, y=PC2, color=tissue))+scale_color_manual(values=custom_palette)+
  theme_bw()+
  xlab("PC3 (11.9%)")+
  ylab("PC2 (16.8%)")+
  std_theme

gempty<-ggplot()+std_theme

pdf(file="../Figures/PCA_123.pdf",
    width=6, height=6, paper="special", bg="white", #used to be 1.75 by 1.75
    fonts="Helvetica", colormodel = "cmyk")
ggarrange(gempty, g21, g31,
          g12, gempty, g32,
          g13, g23, gempty
          , ncol=3, nrow=3)
dev.off()
