library(data.table)
library(tidyverse)
library(patchwork)
library(Biostrings)
library(reshape2)
library(stringr)
library(patchwork)
clonal_families<-c(4,6,8,10,12,14,16,18, 20)
#tools<-c("MiXCR", "changeo","scoper_ID","scoper_hierarchical","scoper_spectral")
SHM<-c("0_001","0_005", "0_01","0_05","0_1","0_2")
leaves <- c("10","20","50","100")
balance <- c("0_0","0_3","0_5","1_0","1_3")
sims<-seq(1,50,1)
path<-"/home1/kavoss/panfs/method_comparison/"

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



for (x in clonal_families) {
  print(x)
  curr_path<-paste(path,as.character(x),sep = "")
  filenames <- list.files(path=curr_path, pattern="analysis_no_singletons.tsv", full.names=TRUE,recursive = TRUE)
  total_df <- do.call(rbind,lapply(filenames,read.csv,sep="\t"))
  
  number_families<-tidyr::expand_grid(unique(total_df$tool),SHM, leaves,balance)
  number_families$clones<-x
  colnames(number_families)<- c("tool","SHM","leaves","balance","clones")
  
  number_families$MSE_num_fam<-apply(number_families,1, FUN=get_MSE_fam_number, df=total_df)
  number_families$MSE_fam_size<-apply(number_families,1, FUN=get_MSE_median_fam_size, df=total_df)
  
  assign(paste("clones",as.character(x),sep = ""),number_families)
  assign(paste("all",as.character(x),sep = ""),total_df)
}


combined<-rbindlist(list(clones4,clones6, clones8,clones10,clones12,clones14,clones16,clones18,clones20))

combined$clones<-as.character(combined$clones)

write.csv(combined, "/home1/kavoss/panfs/method_comparison/output_without_singletons.csv", row.names=FALSE)

all_combined<-rbindlist(list(all4,all6, all8,all10,all12,all14,all16,all18,all20))

write.csv(all_combined, "/home1/kavoss/panfs/method_comparison/output_all_without_singletons.csv", row.names=FALSE)

MSE_SHM_num_fam<-ggplot(combined, aes(SHM,MSE_num_fam, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Number of Families (log)")
MSE_SHM_num_fam


MSE_clones_num_fam<-ggplot(combined, aes(clones,MSE_num_fam, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$clones))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Number of Families (log)")
MSE_clones_num_fam


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


MSE_SHM_size_fam<-ggplot(combined, aes(SHM,MSE_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Median Family Size (log)")
MSE_SHM_size_fam


MSE_clones_size_fam<-ggplot(combined, aes(clones,MSE_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$clones))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Median Family Size (log)")
MSE_clones_size_fam

MSE_leaves_size_fam<-ggplot(combined, aes(leaves,MSE_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$leaves))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Median Family Size (log)")
MSE_leaves_size_fam


MSE_size_fam<-ggplot(combined, aes(tool,MSE_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Median Family Size (log)")
MSE_size_fam

MSE_SHM_num_fam+MSE_clones_num_fam+ MSE_leaves_num_fam+plot_annotation('number of discerned families',theme=theme(plot.title=element_text(hjust=0.5)))+ plot_layout(guides = "collect")
MSE_SHM_size_fam+MSE_clones_size_fam+ MSE_leaves_size_fam+plot_annotation('median size of discerned families',theme=theme(plot.title=element_text(hjust=0.5)))+ plot_layout(guides = "collect")

