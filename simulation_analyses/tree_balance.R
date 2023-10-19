#!/usr/bin/env Rscript
#library("optparse")
library(data.table)
library(tidyverse)
library(apTreeshape)
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/CollessLike/CollessLike_1.0.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
library(CollessLike)

path<-"/home1/kavoss/panfs/method_comparison/"
#folders <- c( "10/0_1/50/", "20/0_001/50/","20/0_1/50/") #"10/0_001/50/"
folders <- c( "10/0_2/50/", "20/0_2/50/") 
balances <- c("0_0","0_3","0_5","1_0","1_3")
sims<-seq(1,50,1)

get_indices<-function(filepath) {
  tree <- read.tree(filepath)
  colless<-colless.like.index(tree, f.size = "ln", diss = "MDM", norm = FALSE)
  sackin <- sackin.index(tree, norm = FALSE)
  return(paste(colless,sackin))
}

total_df<-data.table()
for (x in folders) {
  print(x)
  for (b in balances){
    curr_path<-paste(path,as.character(x),b,sep = "")
    filenames <- list.files(path=curr_path, pattern=".raxml.bestTree$", full.names=TRUE,recursive = TRUE)
    balance_df<-data.frame(filenames)
    balance_df$indices <- lapply(balance_df$filenames,get_indices)
    balance_df<-balance_df %>%
      separate(indices, c("colless", "sackin"), " ")
    balance_df$params<-gsub(path,"",balance_df$filenames)
    balance_df$params<-gsub(x,"",balance_df$params)
    balance_df$balance<- substr(balance_df$params, 0, 3)
    #overview<-tidyr::expand_grid(balances,sims)
    total_df<-rbind(total_df, balance_df)
  }
  
}

total_df$colless<-as.numeric(total_df$colless)
total_df$sackin<-as.numeric(total_df$sackin)

ggplot(total_df, aes(balance,colless)) + 
  scale_y_continuous(trans='log10')+
  geom_boxplot()

ggplot(total_df, aes(balance,sackin)) +
  scale_y_continuous(trans='log10')+
  geom_boxplot()


ggplot(total_df, aes(balance,colless)) + 
  geom_point()

ggplot(total_df, aes(balance,sackin)) + 
  geom_point()



