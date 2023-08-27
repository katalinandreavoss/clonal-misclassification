library(data.table)
library(tidyverse)
library(patchwork)
library(Biostrings)
library(reshape2)
library(stringr)
library(patchwork)
clonal_families<-c(4,6,8,10,12,14,16,18, 20)
tools<-c("SCOPer identicalClones","SCOPer hierarchicalClones", "SCOPer spectralClones")
path<-"/Users/kavoss/Documents/Research/simulations/"


get_values_scoper<-function(dpath,file_name) {
  scoper<-read.table(paste0(dpath,file_name),sep="\t",header=TRUE, fill=TRUE, row.names=NULL)
  scoper_sum<-scoper %>% group_by(clone_id) %>% summarise(n = n())
  number_fams <- length(scoper_sum[scoper_sum$n!=1,]$n)
  med_fam_size <- median(scoper_sum[scoper_sum$n!=1,]$n)
  return(paste(number_fams,med_fam_size))
}

get_MSE_median_fam_size<-function(data,df) {
  clone=data[['clones']]
  shm=data[['SHM']]
  leaf=data[['leaves']]
  which_tool=data[['tool']]
  real<-df%>%
    filter(clones == clone)%>%
    filter(SHM == shm)%>%
    filter(leaves == leaf)%>%
    filter(tool == which_tool)%>%
    select(median_family_size_real)
  pred<-df%>%
    filter(clones == clone)%>%
    filter(SHM == shm)%>%
    filter(leaves == leaf)%>%
    filter(tool == which_tool)%>%
    select(median_family_size)
  mse<-mean((as.numeric(real$median_family_size_real) - as.numeric(pred$median_family_size))^2)
  return(mse)
}

get_MSE_fam_number<-function(data,df) {
  clone=data[['clones']]
  shm=data[['SHM']]
  leaf=data[['leaves']]
  which_tool=data[['tool']]
  real<-df%>%
    filter(clones == clone)%>%
    filter(SHM == shm)%>%
    filter(leaves == leaf)%>%
    filter(tool == which_tool)%>%
    select(real_fam_number)
  pred<-df%>%
    filter(clones == clone)%>%
    filter(SHM == shm)%>%
    filter(leaves == leaf)%>%
    filter(tool == which_tool)%>%
    select(number_families)
  mse<-mean((as.numeric(real$real_fam_number) - as.numeric(pred$number_families))^2)
  return(mse)
}


for (x in clonal_families) {
  print(x)
  curr_path<-paste(path,as.character(x),sep = "")
  filenames <- list.files(path=curr_path, pattern="results_db_idClones.tsv", full.names=TRUE,recursive = TRUE)
  
  total_df<-data.frame(filenames)
  total_df$params<-gsub(path,"",total_df$filenames)
  total_df$params<-gsub("results_db_idClones.tsv","",total_df$params)
  total_df$params<-str_sub(total_df$params,0,-2)
  
  total_df<-total_df %>%
    separate(params, c("clones","SHM","leaves","sim"), "/")
  
  total_df$filenames<-gsub("results_db_idClones.tsv","",total_df$filenames)
  total_df$real_fam_number<-as.numeric(lapply(total_df$filenames, get_real_family_number))
  
  real_freq <-
    filenames <- list.files(path=curr_path, pattern="family_sizes.txt", full.names=TRUE,recursive = TRUE)%>% 
    map_df(~fread(.))
  real_freq$V1<-gsub("_",".",real_freq$V1)
  colnames(real_freq)<-c("clones","SHM","leaves","sim","freq")
  total_df$median_family_size_real<-apply(total_df,1,FUN=get_median_fam_size,df=real_freq)
  
  scoper1_df<-total_df[c("filenames","clones","SHM","leaves","sim","real_fam_number","median_family_size_real")]
  scoper2_df<-total_df[c("filenames","clones","SHM","leaves","sim","real_fam_number","median_family_size_real")]
  scoper3_df<-total_df[c("filenames","clones","SHM","leaves","sim","real_fam_number","median_family_size_real")]
  
  scoper1_df$mixcr<-lapply(scoper1_df$filenames,get_values_scoper,file_name="results_db_idClones.tsv")
  scoper1_df<-scoper1_df %>%
    separate(mixcr, c("number_families", "median_family_size"), " ")
  
  scoper2_df$mixcr<-lapply(scoper2_df$filenames,get_values_scoper,file_name="results_hierClones.tsv")
  scoper2_df<-scoper2_df %>%
    separate(mixcr, c("number_families", "median_family_size"), " ")
  
  scoper3_df$mixcr<-lapply(scoper3_df$filenames,get_values_scoper,file_name="results_specClones.tsv")
  scoper3_df<-scoper3_df %>%
    separate(mixcr, c("number_families", "median_family_size"), " ")
  
  number_families<-tidyr::expand_grid(tools,unique(total_df$clones),unique(total_df$SHM), unique(total_df$leaves))
  colnames(number_families)<- c("tool","clones","SHM","leaves")
  new<-rbindlist(list(scoper1_df,scoper2_df,scoper3_df), idcol = "tool",use.names=TRUE)
  new$tool<-as.character(new$tool)
  new[tool==1]$tool<-"SCOPer identicalClones"
  new[tool==2]$tool<-"SCOPer hierarchicalClones"
  new[tool==3]$tool<-"SCOPer spectralClones"
  
  number_families$MSE_fam_size<-apply(number_families,1, FUN=get_MSE_median_fam_size, df=new)
  
  number_families$MSE_num_fam<-apply(number_families,1, FUN=get_MSE_fam_number, df=new)
  
  #number_families$mixcr.MSE_fam_size<-apply(number_families,1, FUN=get_MSE_median_fam_size, df=total_df, tool="mixcr")
  #number_families$mixcr.MSE_num_fam<-apply(number_families,1, FUN=get_MSE_fam_number, df=total_df, tool="mixcr")
  
  assign(paste("clones",as.character(x),sep = ""),number_families)
  assign(paste("all",as.character(x),sep = ""),new)
  
}

combined<-rbindlist(list(clones4,clones6, clones8,clones10,clones12,clones14,clones16,clones18,clones20))


old<-read.table("/Users/kavoss/Documents/Research/simulations/output.csv",sep=",",header=TRUE, fill=TRUE, row.names=NULL)

combined<-rbindlist(list(old,combined))
combined$clones<-as.character(combined$clones)


MSE_SHM_num_fam<-ggplot(combined, aes(SHM,MSE_num_fam, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Number of Families (log)")

MSE_clones_num_fam<-ggplot(combined, aes(clones,MSE_num_fam, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$clones))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Number of Families (log)")



MSE_leaves_num_fam<-ggplot(combined, aes(leaves,MSE_num_fam, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$leaves))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Number of Families (log)")


MSE_SHM_size_fam<-ggplot(combined, aes(SHM,MSE_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Median Family Size (log)")

MSE_clones_size_fam<-ggplot(combined, aes(clones,MSE_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$clones))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Median Family Size (log)")

MSE_leaves_size_fam<-ggplot(combined, aes(leaves,MSE_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$leaves))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Median Family Size (log)")


MSE_SHM_num_fam+MSE_clones_num_fam+ MSE_leaves_num_fam+plot_annotation('number of discerned families',theme=theme(plot.title=element_text(hjust=0.5)))+ plot_layout(guides = "collect")
MSE_SHM_size_fam+MSE_clones_size_fam+ MSE_leaves_size_fam+plot_annotation('median size of discerned families',theme=theme(plot.title=element_text(hjust=0.5)))+ plot_layout(guides = "collect")
ggsave("/Users/kavoss/Documents/Research/family_size.jpg")
MSE_SHM_num_fam
MSE_clones_num_fam
MSE_leaves_num_fam
