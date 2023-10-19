library('EnvStats')
library(tidyverse)
library(reshape2)
library(dplyr)
gg_qq_empirical <- function(a, b, quantiles = seq(0, 1, 0.01))
{
  a_lab <- deparse(substitute(a))
  if(missing(b)) {
    b <- rnorm(length(a), mean(a), sd(a))
    b_lab <- "normal distribution"
  }
  else b_lab <- deparse(substitute(b))
  
  ggplot(mapping = aes(x = quantile(a, quantiles), 
                       y = quantile(b, quantiles))) + 
    geom_point() +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
    labs(x = paste(deparse(substitute(a)), "quantiles"), 
         y = paste(deparse(substitute(b)), "quantiles"),
         title = paste("Empirical qq plot of", a_lab, "against", b_lab))
}

new_combined<-read.csv("/home1/kavoss/panfs/method_comparison/output_melted.csv")
qqPlot(real_sizes$X18,real_sizes$X18)

path<-"/panfs/qcb-panasas/kavoss/method_comparison/20/0_005/10/0_0/10/"
real_sizes<-read.csv(paste0(path,"/family_sizes.txt"))
colnames(real_sizes)<-c("clones","SHM","leaves","sim","size")

scoper_sizes_hier<-read.csv(paste0(path,"results_hierClones.tsv"),sep="\t")
scoper_sum_hier<-scoper_sizes_hier %>% group_by(clone_id) %>% summarise(n = n())

scoper_sizes_id<-read.csv(paste0(path,"results_db_idClones.tsv"),sep="\t")
scoper_sum_id<-scoper_sizes_id %>% group_by(clone_id) %>% summarise(n = n())

scoper_sizes_sp<-read.csv(paste0(path,"results_specClones.tsv"),sep="\t")
scoper_sum_sp<-scoper_sizes_sp %>% group_by(clone_id) %>% summarise(n = n())

mixcr<-read.csv(paste0(path,"clean.fasta.vdjca.clns_IGH.tsv"),sep="\t")

changeo<-read.csv(paste0(path,"vquest_files/combined_db-pass_clone-pass.tsv"),sep="\t")
changeo_sum<-changeo %>% group_by(clone_id) %>% summarise(n = n())


real<-data.frame(sizes=real_sizes$size)
real$tool<-"real"
scoper_hier<-data.frame(sizes=scoper_sum_hier$n)
scoper_hier$tool<-"scoper_hier"
scoper_id<-data.frame(sizes=scoper_sum_id$n)
scoper_id$tool<-"scoper_id"
scoper_sp<-data.frame(sizes=scoper_sum_sp$n)
scoper_sp$tool<-"scoper_sp"
mixcr_result<-data.frame(sizes=mixcr$readCount)
mixcr_result$tool<-"mixcr_result"
changeo_result<-data.frame(sizes=changeo_sum$n)
changeo_result$tool<-"changeo"



data<-rbind(real,scoper_hier,scoper_sp,scoper_id,mixcr_result,changeo_result)
data<-rbind(real,scoper_hier)

ggplot(data) +geom_qq(aes(sample=sizes, color = tool)) 
ggplot(data[data$sizes!=1,]) +geom_qq(aes(sample=sizes, color = tool)) 

qq <- gg_qq_empirical(real_sizes$X7, scoper_sum_hier$n)
qq + theme_light() + coord_equal()



