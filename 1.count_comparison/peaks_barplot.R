library(tidyr)
library(reshape2)
library(ggplot2)
library(ggthemes)
setwd("~/Oxford 2.0/Scripts/CNN_project/Data")
filter_comp <- read.csv("filter_comparison.csv")
colnames(filter_comp)[c(1,2)] <- c("Stage", "1CPM")
filter_comp <- filter_comp[c("Stage", "Original", "Conservative", "1CPM")]
filter_comp$Stage <- factor(filter_comp$Stage)
filter_comp <- filter_comp[c(8,2,5,7,6,4,3,1), ]
filter_comp$Stage <- factor(filter_comp$Stage, levels=c("iPSC","DE","GT","PF","PE","EP","EN","BLC"))
filter_comp <- melt(filter_comp, id.vars = 1)

ggplot(data=filter_comp, aes(x=Stage, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) + theme_economist_white() + ggtitle("Peaks in .bed-files Before Preprocessing")
  