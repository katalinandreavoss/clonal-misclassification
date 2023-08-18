#!/usr/bin/env Rscript
#library("optparse")
library(data.table)
library(tidyverse)
library(patchwork)
library(Biostrings)
library(reshape2)
library(stringr)
library(patchwork)
clonal_families<-c(4,6,8,10,12,14,16,18, 20)
tools<-c("PTP","MiXCR")
path<-"/Users/kavoss/Documents/Research/simulations/"

get_values_mixcr<-function(filepath) {
  mixcr<-read.table(paste0(filepath,"clean.fasta.vdjca.clns_IGH.tsv"),header=TRUE, fill=TRUE)
  mixcr_sum<-mixcr %>% group_by(readCount) %>% summarise(n = n())
  number_fams <- length(mixcr[mixcr$readCount!=1,]$readCount)
  med_fam_size <- median(mixcr_sum[mixcr_sum$n!=1,]$n)
  return(paste(number_fams,med_fam_size))
}


get_num_fam_no_singletons <- function(file_path) {
  df<-fread(file_path)
  num_fam=df[1,]$V2
  return(num_fam)
}
get_freq_no_singletons <- function(file_path) {
  df<-fread(file_path)
  num_fam=df[2,]$V2
  return(num_fam)
}

get_real_family_number<-function(path) {
  clean_data = readDNAStringSet(paste0(path,"clean.fasta"))
  seq_name = names(clean_data)
  sequence = paste(clean_data)
  clean_data<-data.table(seq_name,sequence)
  clean_data<-data.table(colsplit(clean_data$seq_name, "clone", c("family", "clone")),sequence)
  size<-length(unique(clean_data$family))
  return(size)
}



get_median_fam_size<-function(data,df) {
  clone=data[['clones']]
  shm=data[['SHM']]
  leaf=data[['leaves']]
  s=data[['sim']]
  med<-df%>%
    filter(clones == clone)%>%
    filter(SHM == shm)%>%
    filter(leaves == leaf)%>%
    filter(sim == s)%>%
    select(freq)
  med_real<-median(med$freq)
  return(med_real)
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
  filenames <- list.files(path=curr_path, pattern="mptp_data.txt", full.names=TRUE,recursive = TRUE)
  
  total_df<-data.frame(filenames)
  total_df$params<-gsub(path,"",total_df$filenames)
  total_df$params<-gsub("mptp_data.txt","",total_df$params)
  total_df$params<-str_sub(total_df$params,0,-2)
  
  total_df<-total_df %>%
    separate(params, c("clones","SHM","leaves","sim"), "/")
  
  total_df$number_families<-as.numeric(lapply(total_df$filenames, get_num_fam_no_singletons))
  total_df$median_family_size<-as.numeric(lapply(total_df$filenames, get_freq_no_singletons))
  
  total_df$filenames<-gsub("mptp_data.txt","",total_df$filenames)
  total_df$real_fam_number<-as.numeric(lapply(total_df$filenames, get_real_family_number))
  
  real_freq <-
    filenames <- list.files(path=curr_path, pattern="family_sizes.txt", full.names=TRUE,recursive = TRUE)%>% 
    map_df(~fread(.))
  real_freq$V1<-gsub("_",".",real_freq$V1)
  colnames(real_freq)<-c("clones","SHM","leaves","sim","freq")
  total_df$median_family_size_real<-apply(total_df,1,FUN=get_median_fam_size,df=real_freq)
  
  mixcr_df<-total_df[c("filenames","clones","SHM","leaves","sim","real_fam_number","median_family_size_real")]
  mixcr_df$mixcr<-lapply(mixcr_df$filenames,get_values_mixcr)
  mixcr_df<-mixcr_df %>%
    separate(mixcr, c("number_families", "median_family_size"), " ")
  
  number_families<-distinct(data.frame(tool=tools,clones=x,SHM=total_df$SHM, leaves=total_df$leaves))
  
  new<-rbindlist(list(total_df,mixcr_df), idcol = "tool",use.names=TRUE)
  new$tool<-as.character(new$tool)
  new[tool==1]$tool<-"PTP"
  new[tool==2]$tool<-"MiXCR"
  
  number_families$MSE_fam_size<-apply(number_families,1, FUN=get_MSE_median_fam_size, df=new)
  
  number_families$MSE_num_fam<-apply(number_families,1, FUN=get_MSE_fam_number, df=new)
  
  #number_families$mixcr.MSE_fam_size<-apply(number_families,1, FUN=get_MSE_median_fam_size, df=total_df, tool="mixcr")
  #number_families$mixcr.MSE_num_fam<-apply(number_families,1, FUN=get_MSE_fam_number, df=total_df, tool="mixcr")
  
  assign(paste("clones",as.character(x),sep = ""),number_families)
  assign(paste("all",as.character(x),sep = ""),new)
}

all_combined<-rbindlist(list(all4,all6, all8,all10,all12,all14,all16,all18,all20))

#combined$MSE_SHM<-as.numeric(lapply(combined$SHM, get_MSE_number_families, column_name=quo(SHM),df=all_combined))
#combined$MSE_clones<-as.numeric(lapply(combined$clones, get_MSE_number_families, column_name=quo(clones),df=all_combined))
#combined$MSE_leaves<-as.numeric(lapply(combined$leaves, get_MSE_number_families, column_name=quo(leaves),df=all_combined))

#combined$MSE_SHM_family_size<-as.numeric(lapply(combined$SHM, get_MSE_family_size, column_name=quo(SHM),df=all_combined))
#combined$MSE_clones_family_size<-as.numeric(lapply(combined$clones, get_MSE_family_size, column_name=quo(clones),df=all_combined))
#combined$MSE_leaves_family_size<-as.numeric(lapply(combined$leaves, get_MSE_family_size, column_name=quo(leaves),df=all_combined))


combined<-rbindlist(list(clones4,clones6, clones8,clones10,clones12,clones14,clones16,clones18,clones20))

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

#MSE_SHM_num_fam+MSE_SHM_size_fam+plot_annotation('SHM',theme=theme(plot.title=element_text(hjust=0.5)))
#MSE_clones_num_fam+MSE_clones_size_fam+plot_annotation('clones',theme=theme(plot.title=element_text(hjust=0.5)))
#MSE_leaves_num_fam+MSE_leaves_size_fam+plot_annotation('clones',theme=theme(plot.title=element_text(hjust=0.5)))

all_combined$number_families<-as.numeric(all_combined$number_families)
all_combined$SHM<-as.numeric(gsub("_",".",all_combined$SHM))
all_combined$clones<-as.numeric(all_combined$clones)
all_combined$leaves<-as.numeric(all_combined$leaves)
ggplot(all_combined, aes(real_fam_number,number_families)) + 
  geom_point(aes(shape = tool, color = SHM))+
  scale_color_gradient(low="blue", high="red")


ggplot(all_combined, aes(real_fam_number,number_families)) + 
  geom_point(aes(shape = tool, color = clones))+
  scale_color_gradient(low="blue", high="red")

ggplot(all_combined, aes(real_fam_number,number_families)) + 
  geom_point(aes(shape = tool, color = leaves))+
  scale_color_gradient(low="blue", high="red")


ggplot(all_combined, aes(real_fam_number,number_families)) + 
  geom_point(aes(color = tool))


