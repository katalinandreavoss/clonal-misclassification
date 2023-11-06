library(patchwork)
library(ggplot2)
clonal_families<-c(16)
SHM<-c("0_001","0_005", "0_01","0_05","0_1","0_2")
leaves <- c("10","20","50","100")
balance <- c("0_0")
sims<-seq(1,50,1)
path<-"/scratch1/kavoss/method_comparison/"

listed<-tidyr::expand_grid(clonal_families,SHM, leaves,balance,sims)
listed$filenames<-paste(path,listed$clonal_families,listed$SHM,listed$leaves,listed$balance,listed$sims,"mutations.tsv", sep="/")

filenames <- listed$filenames
total_df <- do.call(rbind,lapply(filenames,read.csv,sep="\t"))

ggplot(total_df, aes(x=SHM,y=mean_mutations))+
  geom_boxplot()
