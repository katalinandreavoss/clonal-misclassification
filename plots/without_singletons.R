library(data.table)
library(tidyverse)
library(patchwork)
library(Biostrings)
library(reshape2)
library(stringr)
library(patchwork)
library(ggplot2)
clonal_families<-c(10,12,14,16,18, 20)
clonal_families<-c(10)
clonal_families<-c(16)
clonal_families<-c(10,20,50)
#tools<-c("MiXCR", "changeo","scoper_ID","scoper_hierarchical","scoper_spectral")
SHM<-c("0_0005", "0_00075","0_001","0_0015","0_002","0_0025","0_003","0_004","0_005", "0_01","0_05","0_1","0_2")
SHM<-c("0_001", "0_005", "0_01","0_05","0_1","0_2")
SHM<-c("0_005", "0_01","0_1","0_2")
SHM<-c("0_01")
leaves <- c("10","20","50","100")
leaves <- c("20","50")
junction_length <- c("10","20","30","40","50","60")
junction_length<-c("20","30","40")
balance<-c("0_0")
sims<-seq(1,10,1)
path<-"/scratch2/kavoss/simulations_junctions"
path<-"/scratch1/kavoss/method_comparison/"
path<-"/panfs/qcb-panasas/kavoss/method_comparison/"
path<-"/scratch1/kavoss/compare_clone_count/"

get_MSE_median_fam_size<-function(data,df) {
  clone=data[['clones']]
  shm=data[['SHM']]
  leaf=data[['leaves']]
  balances=data[['balance']]
  which_tool=data[['tool']]
  real<-df%>%
    filter(clones == clone)%>%
    filter(SHM == shm)%>%
    filter(leaves == leaf)%>%
    filter(balance == balances)%>%
    filter(tool == which_tool)%>%
    select(real_med_size)
  pred<-df%>%
    filter(clones == clone)%>%
    filter(SHM == shm)%>%
    filter(leaves == leaf)%>%
    filter(balance == balances)%>%
    filter(tool == which_tool)%>%
    select(median_size_families)
  mse<-mean((as.numeric(real$real_med_size) - as.numeric(pred[pred$median_size_families!="NaN",]))^2)
  return(mse)
}

get_MSE_fam_number<-function(data,df) {
  clone=data[['clones']]
  shm=data[['SHM']]
  leaf=data[['leaves']]
  balances=data[['balance']]
  which_tool=data[['tool']]
  real<-df%>%
    filter(clones == clone)%>%
    filter(SHM == shm)%>%
    filter(leaves == leaf)%>%
    filter(balance == balances)%>%
    filter(tool == which_tool)%>%
    select(real_num_fam)
  pred<-df%>%
    filter(clones == clone)%>%
    filter(SHM == shm)%>%
    filter(leaves == leaf)%>%
    filter(balance == balances)%>%
    filter(tool == which_tool)%>%
    select(num_families)
  mse<-mean((as.numeric(real$real_num_fam) - as.numeric(pred$num_families))^2)
  return(mse)
}


listed<-tidyr::expand_grid(clonal_families,SHM, leaves,junction_length,sims)
listed$filenames<-paste(path,listed$clonal_families,listed$SHM,listed$leaves,listed$junction_length,listed$sims,"analysis_no_singletons.tsv", sep="/")


filenames <- listed$filenames
total_df <- do.call(rbind,lapply(filenames,read.csv,sep="\t"))

write.csv(total_df, paste0(path,"output_all_without_singletons.csv"), row.names=FALSE)
#total_df<-fread(paste0(path,"16_output_all_without_singletons.csv"))

wrong<-total_df[total_df$tool!="scoper_id",]
wrong_num<-wrong[wrong$num_families!=wrong$real_num_fam,]
wrong_med<-wrong[wrong$median_size_families!=wrong$real_med_size,]
wrong$leaves<-as.factor(wrong$leaves)
wrong$balance<-as.factor(wrong$balance)


ggplot(wrong, aes(SHM,num_families, fill=tool)) + 
  geom_boxplot()+
  ylab("Number of Families")+
  scale_fill_brewer(palette = "Paired")

ggplot(wrong[wrong$SHM=="0_01",], aes(leaves,num_families, fill=tool)) + 
  geom_boxplot()+
  ylab("Number of Families at SHM 0.01")+
  scale_fill_brewer(palette = "Paired")

ggplot(wrong[wrong$leaves==100,], aes(SHM,num_families, fill=tool)) + 
  geom_boxplot()+
  ylab("Number of Families at 100 leaves")+
  scale_fill_brewer(palette = "Paired")

ggplot(wrong, aes(balance,num_families, fill=tool)) + 
  geom_boxplot()+
  ylab("Number of Families")+
  scale_fill_brewer(palette = "Paired")+
  xlab("junction length")

subset<-c("0_0005","0_001","0_0025","0_005", "0_01","0_05","0_1","0_2")

subset_df<-wrong[wrong$SHM %in% subset,]
subset_df$leaves <- as.factor(subset_df$leaves)

ggplot(subset_df[subset_df$tool=="scoper_spec",], aes(leaves,num_families, color=SHM)) + 
  geom_boxplot()+
  ylab("Number of Families")


number_families<-tidyr::expand_grid(unique(total_df$tool),clonal_families,SHM, leaves,junction_length)
colnames(number_families)<- c("tool","clones","SHM","leaves","balance")

number_families$MSE_num_fam<-apply(number_families,1, FUN=get_MSE_fam_number, df=total_df)
number_families$MSE_fam_size<-apply(number_families,1, FUN=get_MSE_median_fam_size, df=total_df)
#number_families$MSE_var_fam_size<-apply(number_families,1, FUN=get_MSE_var_fam_size, df=total_df)
number_families$clones<-as.character(number_families$clones)
number_families$tool <- factor(number_families$tool, levels = c("mixcr", "changeo", "scoper_hier","scoper_spec","mptp"))


combined<-number_families[number_families$tool!="scoper_id",]
#combined<-rbindlist(list(clones4,clones6, clones8,clones10,clones12,clones14,clones16,clones18,clones20))

combined$clones<-as.character(combined$clones)

write.csv(combined, paste0(path,"combined_without_singletons.csv"), row.names=FALSE)

combined<-fread(paste0(path,"combined_without_singletons_scopers.csv"))
#all_combined<-rbindlist(list(all4,all6, all8,all10,all12,all14,all16,all18,all20))

#write.csv(combined, "/scratch1/kavoss/method_comparison/16_output_all_without_singletons.csv", row.names=FALSE)

MSE_SHM_num_fam<-ggplot(combined, aes(SHM,MSE_num_fam, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Number of Families (log)")+
  scale_fill_brewer(palette = "Paired")
MSE_SHM_num_fam


MSE_clones_num_fam<-ggplot(combined, aes(clones,MSE_num_fam, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$clones))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Number of Families (log)")+
  scale_fill_brewer(palette = "Paired")
MSE_clones_num_fam

combined$tool <- factor(combined$tool, levels = c("mixcr", "changeo", "scoper_hier","scoper_spec","mptp"))

MSE_leaves_num_fam<-ggplot(combined, aes(leaves,MSE_num_fam, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$leaves))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Number of Families (log)")
MSE_leaves_num_fam


MSE_balances_num_fam<-ggplot(combined, aes(balance,MSE_num_fam, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$balance))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Number of Families (log)")
MSE_balances_num_fam

combined$tool <- factor(combined$tool, levels = c("mixcr", "changeo", "scoper_hier","scoper_spec","mptp"))

MSE_SHM_size_fam<-ggplot(combined) + 
  geom_boxplot(aes(SHM,MSE_fam_size, fill=tool))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Median Family Size (log)")+
  scale_fill_brewer(palette = "Paired")
one<-MSE_SHM_size_fam+
  theme(text = element_text(face = "bold",size = 20))
one
combined$tool <- factor(combined$tool)

scoper_only<-combined[combined$tool %in% c("gmyc", "gmyc_extend"),]

MSE_SHM_size_fam<-ggplot(scoper_only) + 
  geom_boxplot(aes(SHM,MSE_fam_size, fill=tool))+
  ylab("Mean Squared Error - Median Family Size")+
  scale_fill_brewer(palette = "Paired")
one<-MSE_SHM_size_fam+
  theme(text = element_text(face = "bold",size = 20))
one



MSE_clones_size_fam<-ggplot(combined, aes(clones,MSE_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$clones))+
  ylab("Mean Squared Error - Median Family Size (log)")
MSE_clones_size_fam+
  theme(text = element_text(face = "bold",size = 30))+
  scale_fill_brewer(palette = "Paired")

MSE_leaves_size_fam<-ggplot(combined, aes(leaves,MSE_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$leaves))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Median Family Size (log)")
MSE_leaves_size_fam+
  theme(text = element_text(face = "bold",size = 20))+
  scale_fill_brewer(palette = "Paired")


combined$leaves<-as.factor(combined$leaves)
MSE_leaves_size_fam<-ggplot(combined[combined$SHM=="0_01",], aes(leaves,MSE_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$leaves))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Median Family Size (log)")
MSE_leaves_size_fam+
  theme(text = element_text(face = "bold",size = 20))+
  scale_fill_brewer(palette = "Paired")



MSE_size_fam<-ggplot(combined, aes(tool,MSE_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Median Family Size (log)")
MSE_size_fam+
  theme(text = element_text(face = "bold",size = 20))+
  scale_fill_brewer(palette = "Paired")

MSE_SHM_num_fam+MSE_clones_num_fam+ MSE_leaves_num_fam+plot_annotation('number of discerned families',theme=theme(plot.title=element_text(hjust=0.5)))+ plot_layout(guides = "collect")
MSE_SHM_size_fam+MSE_clones_size_fam+ MSE_leaves_size_fam+plot_annotation('median size of discerned families',theme=theme(plot.title=element_text(hjust=0.5)))+ plot_layout(guides = "collect")

