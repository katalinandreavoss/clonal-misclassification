#install.packages("RRphylo")
#install.packages("ggplot2")
#install.packages("dataframe")
#install.packages("LaplacesDemon")
#install.packages("multimode")
#install.packages("performance")
#install.packages("diptest")
#install.packages("MASS") 
#install.packages("reshape2") 
#install.packages("reshape") 

library(RRphylo)
library(ggplot2)
require(ape)
library(LaplacesDemon)
library(multimode)
library(performance)
library(diptest)
#library(MASS) 
library(reshape2) 



test_tree<-read.tree(file = "/home1/kavoss/simulations/20/0_1/20/1/tree_files/mega_tree_.raxml.bestTree")

matrix<-cophenetic(test_tree)

test_frame<-data.frame(matrix)

family<-subset(test_frame, select=c('family_10_clone_1'))
family$fam<-"other"
family[grep("family_10_", rownames(family)),]$fam<-"same"


#ggplot() + aes(test_frame$time)+ geom_histogram(bins=114, colour=color)
ggplot(data=family,aes(x=family_10_clone_1,fill=fam))+geom_histogram()


#is.unimodal(test_frame[grep("family_2", test_frame$clone),]$time)
dip.test(family[family$fam == "same",]$family_10_clone_1)

#check_multimodal(test_frame[grep("family_2", test_frame$clone),]$time)

values<-c(matrix)
melt_data <- melt(matrix, id = colnames(matrix))
melt_data<-melt_data[melt_data$Var1 != melt_data$Var2, ]   
melt_data$Var2<-gsub("_clone_\\d", "", melt_data$Var2)
melt_data$Var1<-gsub("_clone_\\d", "", melt_data$Var1)
melt_data$family<-"other"
melt_data[melt_data$Var1 == melt_data$Var2, ]$family<-"same"

ggplot(data=melt_data,aes(x=value,fill=family))+geom_histogram()

ggplot() + geom_histogram(aes(values))
