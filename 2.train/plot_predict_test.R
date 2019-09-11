library(precrec)
library(ggplot2)
library(grid)
library(tidyr)
library(plyr)
library(dplyr)

setwd("~/Oxford 2.0/Scripts/CNN_project/Data/better_train/iter1")

act = read.table("final_test_set_act.txt")
features = read.table("8col1CPM_samples.txt")
features = as.character(features$V1)
predictions = read.table("iter1.test.txt")
colnames(predictions) <-
  c("BLC", "DE", "EN", "EP", "GT", "PE", "PF", "iPSC")



score_test <-
  join_scores(
    predictions$BLC,
    predictions$DE,
    predictions$EN,
    predictions$EP
    label_test <- join_labels(act$BLC, act$DE)
    mmdat1 <-
      mmdata(
        score_test,
        label_test,
        modnames = c("BLC", "DE"),
        dsids = c(1, 2)
      )
    
    mmcurves <- evalmod(mmdat1)
    autoplot(mmcurves)
    