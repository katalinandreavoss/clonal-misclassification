library(patchwork)
library(ggplot2)
clonal_families<-c(16)
SHM<-c("0_005", "0_01","0_05","0_1","0_2")
leaves <- c("10","20","50","100")
balance=c("6","12","18")
#balance <- c("0_0")
sims<-seq(1,10,1)
path<-"/scratch1/kavoss/sims_fake_V/"

listed<-tidyr::expand_grid(clonal_families,SHM, leaves,balance,sims)
listed$filenames<-paste(path,listed$clonal_families,listed$SHM,listed$leaves,listed$balance,listed$sims,"sensitivity_precision.tsv", sep="/")

filenames <- listed$filenames
sp_df <- do.call(rbind,lapply(filenames,read.csv,sep="\t"))

write.csv(sp_df, paste0(path,"f1.csv"), row.names=FALSE)
sp_df<-fread(paste0(path,"f1.csv"))
sp_df$tool <- factor(sp_df$tool, levels = c("mixcr", "changeo", "scoper_hier","scoper_sp","gmyc", "gmyc_extend","mptp"))

sp_df$tool <- factor(sp_df$tool, levels = c("mixcr", "changeo", "scoper_hier","scoper_sp","mptp"))

sp_df$tool <- factor(sp_df$tool, levels = c("scoper_sp", "scoper_sp_vj", "scoper_hier"))
sp_df<- sp_df[sp_df$tool %in% c("scoper_hier","scoper_sp"),]

ggplot(sp_df, aes(x=SHM, y=sensitivity, fill = tool))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Paired")

ggplot(sp_df, aes(x=SHM, y=precision, fill = tool))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Paired")
sp_df$leaves<-factor(sp_df$leaves)

two<-ggplot(sp_df, aes(x=SHM, y=f1, fill = tool))+
  geom_boxplot()+
  ylab('F1-Score')+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20))
two

sp_df$clones<-as.character(sp_df$clones)
ggplot(sp_df, aes(x=clones, y=f1, fill = tool))+
  geom_boxplot()+
  ylab('F1-Score')+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20))

sp_df$leaves<-as.character(sp_df$leaves)
ggplot(sp_df, aes(x=leaves, y=FN, fill = tool))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20))

ggplot(sp_df, aes(x=leaves, y=FP, fill = tool))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20))
sp_df$junction_length<-factor(sp_df$junction_length)

sp_df$junction_length<-as.character(sp_df$junction_length)
ggplot(sp_df, aes(x=junction_length, y=f1, fill = tool))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Paired")+
  ylab('F1-Score')+
  theme(text = element_text(face = "bold",size = 20))


sp_df$leaves<-as.factor(sp_df$leaves)

ggplot(sp_df, aes(x=leaves, y=f1, fill = tool))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20))


breaks <- c(-Inf, 60,70,80, Inf)  # Intervals: (-Inf, 60], (60, 90], (90, Inf)

# Create labels for categories
labels <- c("<60", "60-70","70-80" ,">80")

# Create a new column with categories
sp_df$length_category <- cut(sp_df$junction_length_scoper, breaks = breaks, labels = labels, include.lowest = TRUE)
ggplot(sp_df, aes(x=length_category, y=f1, fill = tool))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20))




scoper_comp<-sp_df[sp_df$tool %in% c("scoper_hier","scoper_sp"),]


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
