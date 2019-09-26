### Compare to transcription ###
library("plyr")
library("tidyr")
library("ggplot2")
library("WGCNA")
library("GO.db")
library("pheatmap")
library("reshape2")

setwd("~/Oxford 2.0/Scripts/CNN_project/Data/better_motifs/")

ann_complete <- read.csv("neg_ann_complete.txt", row.names = 1)
RNA_seq <- read.table("31.01.2017.Differentiation_v2.gene.tpm.tsv", header = T)

ann_complete$motifname <- as.character(ann_complete$motifname)
ann_complete$motifname <- sapply(strsplit(ann_complete$motifname, '_'), `[[`, 1)
ann_complete$motifname[!(ann_complete$motifname %in% RNA_seq$GeneName)]

detect_motif_expr <- RNA_seq[RNA_seq$GeneName %in% ann_complete$motifname, ]
detect_motif_expr2 <- RNA_seq[grep("NKX6", RNA_seq$GeneName),]
detect_motif_expr <- rbind(detect_motif_expr, detect_motif_expr2)

#Are they correlated?
cor_matrix_select = matrix()
all_select <- WGCNA::cor(, t(ann_allun), method = "spearman")