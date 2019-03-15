# Convienience function to calculate blast hits after GRAMMY adjustment
library(data.table, quietly = TRUE)
args = commandArgs(trailingOnly = TRUE)
filename=args[1]
df<-fread(filename)
total_blast_hits=as.character(sum(df$AdjustedBlast))
cat(total_blast_hits)
