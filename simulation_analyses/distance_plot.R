library(data.table)
library(tidyverse)
library(patchwork)
library(Biostrings)
library(reshape2)
library(stringr)
library(patchwork)
library(ggplot2)

clonal_families<-c(16)
SHM<-c("0_001","0_005", "0_01","0_05","0_1","0_2")
leaves <- c("10","20","50","100")
#junction_length <- c("10","20","30","40","50","60")
#junction_length <- c("6","12","18","24","30")
junction_length <- c("0_0")
sims<-seq(1,10,1)
path<-"/scratch1/kavoss/method_comparison/"


listed<-tidyr::expand_grid(clonal_families,SHM, leaves,junction_length,sims)
listed$filenames<-paste(path,listed$clonal_families,listed$SHM,listed$leaves,listed$junction_length,listed$sims,"distance.tsv", sep="/")


filenames <- listed$filenames
dist_df <- do.call(rbind,lapply(filenames,read.csv,sep="\t"))

write.csv(dist_df, paste0(path,"/distance.csv"), row.names=FALSE)


ggplot(dist_df, aes(SHM, distance)) +
  geom_boxplot()+
  scale_fill_brewer(palette = "Paired")+
  geom_hline(yintercept=5, linetype="dashed")

