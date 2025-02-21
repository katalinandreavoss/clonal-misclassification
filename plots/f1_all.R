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
balance<-c("20","30","40")
sims<-seq(1,10,1)
path<-"/scratch1/kavoss/restricted_genes/"

clonal_families<-c(16)
SHM<-c("0_05","0_1","0_2","0_3")
leaves <- c("20","50")
balance<-c("0_0")
sims<-seq(1,10,1)
path<-"/scratch1/kavoss/method_comparison/"

listed<-tidyr::expand_grid(clonal_families,SHM, leaves,balance,sims)
listed$filenames<-paste(path,listed$clonal_families,listed$SHM,listed$leaves,listed$balance,listed$sims,"f1_all.tsv", sep="/")

filenames <- listed$filenames
sp_df <- do.call(rbind,lapply(filenames,read.csv,sep="\t"))

write.csv(sp_df, paste0(path,"f1_all.csv"), row.names=FALSE)
sp_df<-fread(paste0(path,"f1_all.csv"))

sp_df$tool <- factor(sp_df$tool, levels = c("mixcr", "changeo", "scoper_hier","scoper_sp","mptp"))
sp_df$leaves<-as.factor(sp_df$leaves)
sp_df$junction_length<-as.character(sp_df$junction_length)


sp_df$f1 <- ifelse(is.na(sp_df$f1), 0, sp_df$f1)

ggplot(sp_df, aes(x=SHM, y=sensitivity, fill = tool))+
  geom_boxplot()+
  ylab('F1-Score')+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20),
        panel.background = element_rect(fill = "white") )

ggplot(sp_df, aes(x=SHM, y=f1, fill = tool))+
  geom_boxplot()+
  ylab('F1-Score')+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20),
        panel.background = element_rect(fill = "white") )

ggplot(sp_df, aes(x=SHM, y=singletons, fill = tool))+
  geom_boxplot()+
  ylab('# singletons')+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20),
        panel.background = element_rect(fill = "white") )

ggplot(sp_df[sp_df$SHM =="0_05",], aes(x=leaves, y=f1, fill = tool))+
  geom_boxplot()+
  ylab('F1-Score')+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20),
        panel.background = element_rect(fill = "white") )

ggplot(sp_df, aes(x=leaves, y=f1, fill = tool))+
  geom_boxplot()+
  ylab('F1-Score')+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20),
        panel.background = element_rect(fill = "white") )


ggplot(sp_df, aes(x=leaves, y=unsolved, fill = tool))+
  geom_boxplot()+
  ylab('# unsolved samples per tool')+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20),
        panel.background = element_rect(fill = "white") )

ggplot(sp_df, aes(x=leaves, y=solved, fill = tool))+
  geom_boxplot()+
  ylab('# solved samples per tool')+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20),
        panel.background = element_rect(fill = "white") )

ggplot(sp_df, aes(x=SHM, y=ratio, fill = tool))+
  geom_boxplot()+
  ylab('ratio of solved samples/all samples log')+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20),
        panel.background = element_rect(fill = "white") )



ggplot(sp_df, aes(x=junction_length, y=ratio, fill = tool))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Paired")+
  ylab('ratio')+
  xlab('# mutations added to known V genes to create new V genes')+
  theme(text = element_text(face = "bold",size = 20),
        panel.background = element_rect(fill = "white"))

ggplot(sp_df, aes(x=junction_length, y=f1, fill = tool))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Paired")+
  ylab('ratio')+
  xlab('# mutations added to known V genes to create new V genes')+
  theme(text = element_text(face = "bold",size = 20),
        panel.background = element_rect(fill = "white"))


sp_df$rand <- ifelse(is.na(sp_df$rand), 0, sp_df$rand)

ggplot(sp_df, aes(x=SHM, y=rand, fill = tool))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Paired")+
  ylab('Rand Index')+
  theme(text = element_text(face = "bold",size = 20),
        panel.background = element_rect(fill = "white"))


