### General description ###
### Compare WGCNA module stage influence results to
### CNN filter stage influence and find overlapping TF motifs
### associated with both. 
### End of description ###

library("plyr")
library("tidyr")
library("ggplot2")
library("WGCNA")
library("GO.db")
library("pheatmap")
library("reshape2")

# Get tables from motif directory
setwd("~/Oxford 2.0/Scripts/CNN_project/Data/better_motifs/")
ann_qselect <-
  read.csv(
    file = "ann_qselect.txt",
    header = T,
    sep = ",",
    row.names = 1
  )
ann_qselect2 <-
  read.csv(
    file = "ann_qselect2.txt",
    header = T,
    sep = ",",
    row.names = 1
  )
ann_qselect3 <-
  read.csv(
    file = "ann_qselect3.txt",
    header = T,
    sep = ",",
    row.names = 1
  )
ann_allun <-
  read.csv(
    file = "all_unannotated.txt",
    header = T,
    sep = ",",
    row.names = 1
  )
filter_q_and_motifs <-
  read.csv(
    file = "filter_q_and_motifs.txt",
    header = T,
    sep = ',',
    row.names = 1
  )

# Get the HOMER specific data from "comparison" directory and work there:
setwd("~/Oxford 2.0/Scripts/CNN_project/Data/better_comparison/")
recoding <- read.table(file = "WGCNA_recoding.txt", header = TRUE)
module_eigen <-
  read.table(file = "module_eigengenes.txt", header = T, sep = "")
motif_matches <-
  read.delim(file = "novo_motif_matches.txt", sep = "\t", header = F)
motif_table <-
  read.table(file = "novo_motif_table.txt", header = T, sep = "")

# Adapt filter_q_and_motifs column for easier matching
filter_q_and_motifs$Query.ID <-
  gsub("ilter", "", x = filter_q_and_motifs$Query.ID)

# Correct ann_allun rownames: remove space
row.names(ann_allun) <- gsub(" ", "", row.names(ann_allun))

# Fix the motif_matches table, because the colnames are inappropriate. You need this for merging later
motif_matches$V4 <- NULL
colnames(motif_matches)[4:13] <-
  sapply(motif_matches[1, 3:12], as.character)
colnames(motif_matches)[1:2] <-
  sapply(motif_matches[1, 1:2], as.character)
motif_matches$V3 <- NULL
colnames(motif_matches)[3] <- "database_rank"


#Take averages of the values per stage from Marta's module eigenvalues
module_eigen$stage <-
  sapply(strsplit(row.names(module_eigen), '-'), function(x)
    x[1])
module_eigen$stage <- as.factor(module_eigen$stage)
agg_mean <-
  aggregate(module_eigen[, 1:13],
            by = list(module_eigen$stage),
            FUN = mean)
row.names(agg_mean) <- agg_mean[, 1]
agg_mean$Group.1 <- NULL
module_agg <- as.data.frame(t(agg_mean))
module_agg <- module_agg[, c(6, 2, 5, 8, 7, 4, 3, 1)]

# Redo, with filter names only (so remove TF names), this means you need a smaller table first:
ann_slimmed <-
  sapply(strsplit(row.names(ann_qselect), ' '), function(x)
    x[1])
ann_qselect$row <- ann_slimmed
ann_slimdf <- ann_qselect[!duplicated(ann_slimmed), ]
row.names(ann_slimdf) <-
  sapply(strsplit(row.names(ann_slimdf), ' '), function(x)
    x[1])
ann_slimdf$row <- NULL
ann_qselect$row <- NULL



#Redo, remove irrelevant modules and save
module_agg_select <-
  module_agg[c("black", "green", "red", "turquoise"), ]
write.csv(module_agg_select, "module_agg_select.txt", row.names = TRUE)

# Determine correlation and plot
cor_matrix_select = matrix()
all_select <-
  WGCNA::cor(t(module_agg_select), t(ann_allun), method = "spearman")

pheatmap(
  all_select,
  color = colorRampPalette(c("navy", "white", "red"))(30),
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  show_colnames =  FALSE
)

###################################################################

### Some stats:
# A list of all significant TF binding motifs (q < 0.05):
test <- sapply(strsplit(row.names(ann_qselect), '_'), `[[`, 1)
# How many are unique after clipping after '_'? Ans: 34
length(unique(test))


# For each row in the resulting matrix (all_select), how many are correlated above 0.8?
#Select these
count_df <- data.frame("module", "correlated>0.8")

count_df <- count_df[FALSE,]
for (i in row.names(all_select)) {
  y <- all_select[i, all_select[i, ] > 0.8]
  len_list <- length(y)
  y <- list(y)
  y <- unlist(y)
  z <- names(unlist(y))
  assign(paste("matchlist", i, sep = ""), y)
  assign(paste("naam", i, sep = ""), z)
  new_row <- cbind(i, len_list)
  print(new_row)
  count_df <- rbind(count_df, new_row)
}
colnames(count_df) <- c("module", "correlated>0.8")


### For overlap comparison: ###
# 1) list of TFs per HOMER module
# 2) list of filters correlated above 0.8 with a HOMER module
# 3) list of TFs associated with respective filters

#Prepare an empty df you can fill with annotated, correlated filters
correlated_filters <- data.frame()

# Once you have list of correlated filters, find how many filters are
# Annotated by Tomtom
module_cor_filtname <- ls()[grep("naam", ls())]
for (i in module_cor_filtname) {
  filt_list <- get(i)[get(i) %in% filter_q_and_motifs$Query.ID]
  i2 <- gsub("naam", "", i)
  assign(paste("ann_filt_list", i2, sep = "_"), filt_list)

  # Use the annotated filters to make a table of
  # correlated AND annotated filters, their associated motifs, q-value and the module
  for (j in filt_list) {
    corann_filt_data <- filter_q_and_motifs[filter_q_and_motifs[, 1] == j, ]
    corann_filt_data$cor_module <- i2
    correlated_filters <-
      rbind(correlated_filters, corann_filt_data)

  }
}

# Add a column with correlation values
# This creates a table of filters that are correlated to a module, the q-value of its annotation,
# the TF motif annotation, and the correlation. Save this "master table".
# Warning! this table still contains rows with TF motifs that are not overlapping with the
# WGCNA module tables.
correlated_filters$correlation <- NA

for (i in (1:nrow(correlated_filters))) {
  correlated_filters[i, "correlation"] <-
    all_select[correlated_filters[i, "cor_module"], correlated_filters[i, "Query.ID"]]
}

write.table(correlated_filters, "correlated_ann_filters_cor.txt")

# Determine overlap [inefficient coding]:
black_select <-
  correlated_filters[correlated_filters$cor_module == "black", ]
black_f_in_m <-
  black_select[black_select$motifname %in% motif_matches[motif_matches$module == "black", "ID"], ]

green_select <-
  correlated_filters[correlated_filters$cor_module == "green", ]
green_f_in_m <-
  green_select[green_select$motifname %in% motif_matches[motif_matches$module == "green", "ID"], ]
unique(green_f_in_m$Query.ID)

red_select <-
  correlated_filters[correlated_filters$cor_module == "red", ]
red_f_in_m <-
  red_select[red_select$motifname %in% motif_matches[motif_matches$module == "red", "ID"], ]


turquoise_select <-
  correlated_filters[correlated_filters$cor_module == "turquoise",]
turquoise_f_in_m <-
  turquoise_select[turquoise_select$motifname %in% motif_matches[motif_matches$module == "turquoise", "ID"],]
unique(turquoise_f_in_m$Query.ID)

# A df with all info. Only contains rows with TF motifs that overlap between filters and WGCNA modules.
mod_CNN_TFoverlap <-
  rbind(black_f_in_m, green_f_in_m, red_f_in_m, turquoise_f_in_m)

write.table(mod_CNN_TFoverlap, "mod_CNN_TFoverlap.txt")

################################################################################\


### Developmental and MODY genes in overlap ###
#Check which MODY related genes are detected in the filter associated TF motifs AND the WGCNA modules.
mody_overlap <-
  mod_CNN_TFoverlap[grep(
    "GCK|GLIS3|HNF1A|HNF1B|HNF4A|KCNJ11|PPARG|WFS1|INS|MNX1|NEUROG3|SLC2A2|NKX2",
    mod_CNN_TFoverlap$motifname
  ), ]


# Also check the same for pancreatic developmental genes.
# Import them from a csv compiled from papers, correct where needed:
dev_TFs <- read.csv("common_dev_motifs.csv", header = F)
dev_TFs <- as.character(dev_TFs$V1)
dev_TFs[1] <- "POU5F1"
dev_TFs[2] <- "OCT4"
dev_TFs[32] <- "NKX6.1"
dev_TFs[44] <- "HNF1B"
# Add genes from MODY list (Slides Marta) and remove duplicates
dev_TFs <-
  c(
    dev_TFs,
    "PDX1",
    "PTF1A",
    "GATA4",
    "GATA6",
    "NEUROG3",
    "GLIS3",
    "PAX6",
    "MNX1",
    "RFX6",
    "NEUROD1",
    "NKX2",
    "HNF1B",
    "SOX9",
    "NEK8",
    "NPHP3",
    "UBR1"
  )
dev_TFs <- unique(dev_TFs)

# Once collapsed all into a single string, check if dev. genes occur in the
# WGCNA module vs. CNN filter overlap df:
dev_TFs <- paste(dev_TFs, collapse = " ")
dev_TFs <- gsub(" ", "|", dev_TFs)
dev_overlap <- mod_CNN_TFoverlap[grep(dev_TFs,
                                      mod_CNN_TFoverlap$motifname),]

# Use the dev_TFs df to show how overlap between WGCNA modules and CNN filters contains expected
# dev genes.
write.table(dev_overlap, "dev_overlap.txt")
