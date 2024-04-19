library(data.table)
library(tidyverse)
library(patchwork)
library(Biostrings)
library(reshape2)
library(stringr)
library(patchwork)
library(ggplot2)
clonal_families<-c(16)

SHM<-c("0_001", "0_005", "0_01","0_05","0_1","0_2")
leaves <- c("10","20","50","100")
balance<-c("0_0")
sims<-seq(1,50,1)

path<-"/scratch1/kavoss/method_comparison/"

listed<-tidyr::expand_grid(clonal_families,SHM, leaves,balance,sims)
listed$filenames<-paste(path,listed$clonal_families,listed$SHM,listed$leaves,listed$balance,listed$sims,"analysis_complete.tsv", sep="/")

filenames <- listed$filenames
total_df <- do.call(rbind,lapply(filenames,read.csv,sep="\t"))


write.csv(total_df, paste0(path,"output_all_with_singletons_gmyc_2.csv"), row.names=FALSE)

total_df<-fread(paste0(path,"output_all_with_singletons_gmyc_2.csv"))



total_df<-total_df[!total_df$tool %in% c("scoper_id","gmyc", "gmyc_extend"),]

total_df$tool <- factor(total_df$tool, levels = c("mixcr", "changeo", "scoper_hier","scoper_spec","mptp"))
total_df$leaves<-as.factor(total_df$leaves)

ggplot(total_df, aes(SHM,num_families, fill=tool)) + 
  geom_boxplot()+
  ylab("Number of Families")+
  scale_y_continuous(trans='log10')+
  scale_fill_brewer(palette = "Paired")+
  geom_hline(yintercept=16,linetype="dashed")+
  theme(text = element_text(face = "bold",size = 20))

ggplot(total_df, aes(leaves,num_families, fill=tool)) + 
  geom_boxplot()+
  ylab("Number of Families")+
  scale_y_continuous(trans='log10')+
  scale_fill_brewer(palette = "Paired")+
  geom_hline(yintercept=16,linetype="dashed")

ggplot(total_df, aes(tool,num_families, fill=tool)) + 
  geom_boxplot()+
  ylab("Number of Families")+
  scale_y_continuous(trans='log10')+
  scale_fill_brewer(palette = "Paired")+
  geom_hline(yintercept=16,linetype="dashed")


total_df_without<-fread(paste0(path,"output_all_without_singletons_gmyc_2.csv"))
total_df_without<-total_df_without[!total_df_without$tool %in% c("scoper_id","gmyc", "gmyc_extend"),]

total_df_without$tool <- factor(total_df_without$tool, levels = c("mixcr", "changeo", "scoper_hier","scoper_spec","mptp"))

ggplot(total_df_without, aes(SHM,num_families, fill=tool)) + 
  geom_boxplot()+
  ylab("Number of Families")+
  scale_y_continuous(trans='log10')+
  scale_fill_brewer(palette = "Paired")+
  geom_hline(yintercept=16,linetype="dashed")

ggplot(total_df_without, aes(tool,num_families, fill=tool)) + 
  geom_boxplot()+
  ylab("Number of Families")+
  scale_y_continuous(trans='log10')+
  scale_fill_brewer(palette = "Paired")+
  geom_hline(yintercept=16,linetype="dashed")

