library(patchwork)
library(ggplot2)
clonal_families<-c(16)
SHM<-c("0_001","0_005", "0_01","0_05","0_1","0_2")
leaves <- c("10","20","50","100")
balance <- c("0_0")
sims<-seq(1,50,1)
path<-"/scratch1/kavoss/method_comparison/"

listed<-tidyr::expand_grid(clonal_families,SHM, leaves,balance,sims)
listed$filenames<-paste(path,listed$clonal_families,listed$SHM,listed$leaves,listed$balance,listed$sims,"seq_similarity.tsv", sep="/")

get_sub_df<-function(filename) {
  df<-read.csv(filename,sep="\t")
  families<-Reduce(intersect, list(df[df$tool=="correct",]$family,df[df$tool=="mixcr",]$family,df[df$tool=="scoper_hier",]$family,
                                   df[df$tool=="scoper_sp",]$family,df[df$tool=="changeo",]$family,df[df$tool=="mptp",]$family))
  return_df<-df[df$family %in% families,] %>% distinct()
  cols<-c("family","tool","hamming_dist","lev_dist","seq_similarity","seq_similarity_lev","midpoint_root")
  return_df<-return_df[cols]
  
  combined_df <- return_df %>%
    group_by(family,tool, midpoint_root) %>%
    summarise(across(c("hamming_dist","lev_dist","seq_similarity","seq_similarity_lev"), mean))
  
  #for(family in families){
  #  if(nrow(combined_df[combined_df$family==family,])!=12){
  #    print(filename)
  #    print(family)
  #  }
  #}
  return(combined_df)
}


filenames <- listed$filenames
total_df <- do.call(rbind,lapply(filenames,get_sub_df))



#ggplot(total_df, aes(x=seq_similarity))+
#  geom_histogram()+
#  aes(fill = as.factor(tool))+
#  facet_grid(tool ~ .)
tools<-c("correct","mixcr","changeo","scoper_hier","scoper_sp","mptp")

#ks.test(total_df[total_df$tool=="correct",]$seq_similarity,total_df[total_df$tool=="mixcr",]$seq_similarity)
write.csv(total_df, "/scratch1/kavoss/method_comparison/ancestral_seq.csv", row.names=FALSE)


plot<-ggplot(transform(total_df[total_df$tool %in% tools,],
                 tool=factor(tool,levels=c("correct","mixcr","changeo","scoper_hier","scoper_sp","mptp"))))+
  geom_histogram(aes(seq_similarity),binwidth = 0.25)+
  aes(fill = as.factor(tool))+
  facet_grid(midpoint_root ~ tool )+
  guides(fill=guide_legend(title="Tool"))+
  theme(text = element_text(size=15)) +
  scale_fill_manual(values=c("black","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C","#FB9A99"))

plot

midpoint_df<-total_df[total_df$midpoint_root =="True",]

ggplot(transform(midpoint_df[midpoint_df$tool %in% tools,],
                       tool=factor(tool,levels=c("correct","mixcr","changeo","scoper_hier","scoper_sp","mptp"))))+
  geom_histogram(aes(seq_similarity),binwidth = 0.25)+
  aes(fill = as.factor(tool))+
  facet_grid(. ~ tool )+
  guides(fill=guide_legend(title="Tool"))+
  theme(text = element_text(size=15)) +
  xlim(90,100)+
  scale_fill_manual(values=c("black","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C","#FB9A99"))



total_counts <- total_df[total_df$tool %in% tools,] %>%
  group_by(tool, midpoint_root) %>%
  summarise(total_count = n())

plot + geom_text(data = total_counts,
                 aes(x = Inf, y = Inf, label = total_count),
                 hjust = 1, vjust = 2, size = 4)+
  xlim(90,100)


ggplot(transform(total_df,
                 tool=factor(tool,levels=c("correct","correct_parsim","correct_ml","mixcr","changeo","scoper_hier","scoper_sp","mptp"))))+
  geom_histogram(aes(seq_similarity_lev))+
  aes(fill = as.factor(tool))+
  facet_grid(tool ~ .)+
  scale_fill_manual(values=c("black","darkgrey","lightgrey","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C","#FB9A99"))

ggplot(transform(total_df,
                 tool=factor(tool,levels=c("correct","correct_parsim","correct_ml","mixcr","changeo","scoper_hier","scoper_sp"))))+
  geom_histogram(aes(lev_dist))+
  aes(fill = as.factor(tool))+
  facet_grid(tool ~ .)+
  scale_fill_manual(values=c("black","lightgrey","darkgrey","#A3A500","#F8766D","#00BF7D","#E76BF3"))

ggplot(transform(total_df,
                 tool=factor(tool,levels=c("correct","correct_parsim","correct_ml","mixcr","changeo","scoper_hier","scoper_sp"))))+
  geom_histogram(aes(seq_similarity))+
  aes(fill = as.factor(tool))+
  facet_grid(tool ~ leaves)+
  scale_fill_manual(values=c("black","lightgrey","darkgrey","#A3A500","#F8766D","#00BF7D","#E76BF3"))


#ggplot(transform(total_df,
#                 tool=factor(tool,levels=c("correct","correct_parsim","correct_ml","mixcr","changeo","scoper_hier","scoper_sp"))))+
#  geom_boxplot(aes(seq_similarity))+
#  aes(fill = as.factor(tool))+
#  facet_grid(tool ~ leaves)+
#  scale_fill_manual(values=c("black","lightgrey","darkgrey","#A3A500","#F8766D","#00BF7D","#E76BF3"))


ggplot(total_df)+geom_point(aes(x=hamming_dist, y=lev_dist))


number_families<-tidyr::expand_grid(unique(total_df$tool),clonal_families,SHM, leaves,balance)
colnames(number_families)<- c("tool","clones","SHM","leaves","balance")





df <- data.frame(x=1:5, y=1, col=letters[1:5])

# Construct the plot
g <- ggplot(df, aes(x=x, y=y, color=col)) + 
  geom_point(size=5)+
  scale_colour_brewer(palette = "Paired")

colors <- ggplot_build(g)$data[[1]]$colour




# look at mutations in 0_001

ggplot(subset(total_df, subset=(tool=="correct"&balance=="0_0")))+
  geom_boxplot(aes(x=SHM,y=hamming_dist))
