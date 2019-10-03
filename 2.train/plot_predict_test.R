### This script makes plots for the test set results ###
# After using Basset_predict.py to predict the classes 
# of the test set sequences. Use the predictions and the 
# Original labels to make PR ad ROC plots.
# Load packages
library(precrec)
library(ggplot2)
library(grid)
library(tidyr)
library(plyr)
library(dplyr)

# Go to the data directory with test set results
setwd("~/Oxford 2.0/Scripts/CNN_project/Data/better_train/iter1")

# Get the tables of labels (acutal classes), predictions and overview of the samples
act = read.table("final_test_set_act.txt")
features = read.table("8col1CPM_samples.txt")
features = as.character(features$V1)
predictions = read.table("iter1.test.txt")
colnames(predictions) <-
  c("BLC", "DE", "EN", "EP", "GT", "PE", "PF", "iPSC")

# Prepare predictions and actual labels, using precrec
score_test <-
  join_scores(
    predictions$iPSC,
    predictions$DE,
    predictions$GT,
    predictions$PF,
    predictions$PE,
    predictions$EP,
    predictions$EN,
    predictions$BLC
  )

label_test <-
  join_labels(act$iPSC, act$DE, act$GT, act$PF, act$PE, act$EP, act$EN, act$BLC)

mmdat1 <-
  mmdata(score_test,
         label_test,
         modnames = c("iPSC","DE","GT","PF","PE","EP","EN","BLC"),
         dsids = c(1, 2, 3, 4, 5, 6, 7, 8))

mmcurves <- evalmod(mmdat1)

# Plot the ROC and PR curves for each of the 8 stages in one plot
autoplot(mmcurves, curvetype = "ROC")
autoplot(mmcurves, curvetype = "PRC")
    


