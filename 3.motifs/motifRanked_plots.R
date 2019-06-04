###   Quick Description   ###
#Use this script to make barplots of influence and also lineplots of influence trhoughout development

library(dplyr)
library(fitdistrplus)
library(logspline)
library(ggplot2)
library(reshape2)
library(qvalue)
library(ape)
library(cluster) 
library(ggdendro)
library(dendextend)
library(psych)
library(gplots)
library(wrapr)

setwd("~/Oxford 2.0/Scripts/CNN_project/Data/motif_original_random")

#Procedure for model, read table and select annotated rows. Select top 20.
table_standard <- read.table("table.txt")
standard_annon <- table_standard[table_standard$annotation != ".", ]
standard_annon <- standard_annon[order(standard_annon$std, decreasing= T), ]
standard_annon$annotation <- make.unique(as.character(standard_annon$annotation), sep = "_")
standard_annon$annotation <- factor(standard_annon$annotation, levels =  standard_annon$annotation[order(-standard_annon$std)])
standard_annon <- standard_annon[1:10,]

#Plot
g <- ggplot(standard_annon, aes(x =annotation, y = std))
g + geom_col(aes(fill=annotation)) + theme_minimal() + ggtitle("Annotated filter influence") +
  xlab("Annotation") + ylab("Standard deviation") + theme(plot.title = element_text(hjust = 0.5))



#Import the influence per stage
infl_standard <- read.table("table_target.txt")
infl_standard <- subset(infl_standard, V1 %in%  rownames(standard_annon))
infl_standard$V3 <- ordered(infl_standard$V3, levels = c("iPSC", "DE", "PGT", "PFG", "PE", "EP", "EN", "BLC"))
infl_standard <- infl_standard %>% mutate(
  V1 = plyr::mapvalues(V1,
                       from = rownames(standard_annon), to = as.character(standard_annon$annotation)))

p<-ggplot(infl_standard, aes(x=V3, y=V4, group=V1)) +
  geom_line(aes(color=V1), size =1)+
  geom_point(aes(color=V1), size =2.5) +
  ylim(-0.03,0.08)
p <- p + 
  labs(title="Filter influence per stage",x="Stage", y = "Influence") +
  theme_minimal(base_size = 22)
p + theme(legend.position="top") + theme(plot.title = element_text(hjust = 0.5))


