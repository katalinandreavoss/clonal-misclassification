library(data.table)
library(tidyverse)
library(patchwork)
library(Biostrings)
library(reshape2)
library(stringr)
library(patchwork)
library(ggplot2)

clonal_families<-c(16)
SHM<-c("0_001","0_005", "0_01","0_05","0_1","0_2")
leaves <- c("10","20","50","100")
#junction_length <- c("10","20","30","40","50","60")
#junction_length <- c("6","12","18","24","30")
junction_length <- c("0_0")
sims<-seq(1,10,1)
path<-"/scratch1/kavoss/method_comparison/"

#sims<-sims[ !sims == 28]

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}



listed<-tidyr::expand_grid(clonal_families,SHM, leaves,junction_length,sims)
listed$filenames<-paste(path,listed$clonal_families,listed$SHM,listed$leaves,listed$junction_length,listed$sims,"ham_distance.tsv", sep="/")


filenames <- listed$filenames
ham_dist_df <- do.call(rbind,lapply(filenames,read.csv,sep="\t"))


subset_df <- subset(sp_df, select = c(SHM, leaves, sim,junction_length, junction_length_scoper))
subset_df<-unique(subset_df)

merged_df <- merge(ham_dist_df, subset_df, by = c("SHM", "leaves", "sim", "junction_length"), all.x = TRUE)
write.csv(ham_dist_df, paste0(path,"/hamming.csv"), row.names=FALSE)


ggplot(ham_dist_df, aes(SHM, ham_dist, fill = category)) +
  geom_split_violin()+
  scale_fill_brewer(palette = "Paired")+
  ylim(c(0,0.3))

#merged_df$junction_length_scoper<-as.factor(merged_df$junction_length_scoper)
breaks <- c(-Inf, 60,70,80, Inf)  # Intervals: (-Inf, 60], (60, 90], (90, Inf)

# Create labels for categories
labels <- c("<60", "60-70","70-80" ,">80")

# Create a new column with categories
merged_df$length_category <- cut(merged_df$junction_length_scoper, breaks = breaks, labels = labels, include.lowest = TRUE)

ggplot(merged_df, aes(length_category, ham_dist, fill = category)) +
  geom_split_violin()+
  scale_fill_brewer(palette = "Paired")+
  ylim(c(0,0.5))+
  geom_hline(yintercept=0.15, linetype="dashed")


ggplot(ham_dist_df, aes(junction_length, ham_dist, fill = category)) +
  geom_boxplot()+
  scale_fill_brewer(palette = "Paired")+
  geom_hline(yintercept=0.15, linetype="dashed")



sub_df<-subset(merged_df, SHM == '0_2')#
ggplot(sub_df, aes(length_category, ham_dist, fill = category)) +
  geom_split_violin()+
  scale_fill_brewer(palette = "Paired")+
  ylim(c(0,0.2))+
  ylab("hamming distance at SHM 0.2 with 100 leaves")+
  geom_hline(yintercept=0.15, linetype="dashed")


