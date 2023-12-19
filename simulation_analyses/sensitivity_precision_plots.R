library(patchwork)
library(ggplot2)
clonal_families<-c(16)
SHM<-c("0_001","0_005", "0_01","0_05","0_1","0_2")
leaves <- c("10","20","50","100")
balance <- c("0_0")
sims<-seq(1,50,1)
path<-"/scratch1/kavoss/method_comparison/"

listed<-tidyr::expand_grid(clonal_families,SHM, leaves,balance,sims)
listed$filenames<-paste(path,listed$clonal_families,listed$SHM,listed$leaves,listed$balance,listed$sims,"sensitivity_precision.tsv", sep="/")

filenames <- listed$filenames
total_df <- do.call(rbind,lapply(filenames,read.csv,sep="\t"))

write.csv(total_df, "/scratch1/kavoss/method_comparison/f1.csv", row.names=FALSE)

total_df$tool <- factor(total_df$tool, levels = c("mixcr", "changeo", "scoper_hier","scoper_sp","mptp","gmyc"))

ggplot(total_df, aes(x=SHM, y=sensitivity, fill = tool))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Paired")

ggplot(total_df, aes(x=SHM, y=precision, fill = tool))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Paired")

two<-ggplot(total_df, aes(x=leaves, y=f1, fill = tool))+
  geom_boxplot()+
  ylab('F1-Score')+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20))
two

ggplot(total_df, aes(x=SHM, y=FN, fill = tool))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20))

total_df$leaves<-as.factor(total_df$leaves)

ggplot(total_df, aes(x=leaves, y=f1, fill = tool))+
  geom_boxplot()+
  scale_fill_manual(values=c("#F8766D","#A3A500","#00BF7D","#E76BF3"))

scoper_comp<-total_df[total_df$tool %in% c("scoper_hier","scoper_sp"),]


ggplot(scoper_comp, aes(x=SHM, y=FN, fill = tool))+
  geom_boxplot()+
  ylab('FN')+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20))
ggplot(scoper_comp, aes(x=SHM, y=FP, fill = tool))+
  geom_boxplot()+
  ylab('FP')+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20))

ggplot(scoper_comp, aes(x=SHM, y=f1, fill = tool))+
  geom_boxplot()+
  ylab('F1-Score')+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20))
