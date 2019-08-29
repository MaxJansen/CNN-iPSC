### General description ###
### Compare homer stage influence results to
### CNN (only annotated with q value < 0.05) filter stage influence 
### Also, only filters with a general influence (std) > 0.1
### End of description ###

library("plyr")
library("tidyr")
library("ggplot2")
library("WGCNA")
library("pheatmap")
library("reshape2")

setwd("~/Oxford 2.0/Scripts/CNN_project/Data/homer_comparison/")

recoding <- read.table(file = "WGCNA_recoding.txt", header = TRUE)
ann_qselect <- read.csv(file ="ann_qselect.txt", header = T, sep = "")
module_eigen <- read.table(file = "module_eigengenes.txt", header = T, sep = "")
motif_matches <- read.delim(file = "novo_motif_matches.txt", sep = "\t", header = F)
motif_table <- read.table(file = "novo_motif_table.txt", header = T, sep = "")

#Fix the motif_matches table, because the colnames are inappropriate. You need this for merging later
motif_matches$V4 <- NULL
colnames(motif_matches)[4:13] <- sapply(motif_matches[1,3:12], as.character)
colnames(motif_matches)[1:2] <- sapply(motif_matches[1,1:2], as.character)
motif_matches$V3 <- NULL
colnames(motif_matches)[3] <- "database_rank"

#Take averages of the values per stage from Marta's module eigenvalues
module_eigen$stage <- sapply(strsplit(row.names(module_eigen), '-'), function(x) x[1])
module_eigen$stage <- as.factor(module_eigen$stage)
agg_mean <- aggregate(module_eigen[,1:13],by=list(module_eigen$stage),FUN=mean)
row.names(agg_mean) <- agg_mean[,1]
agg_mean$Group.1 <- NULL
module_agg <- as.data.frame(t(agg_mean))
module_agg<- module_agg[,c(6,2,5,8,7,4,3,1)]

# Calculate the correlations and plot
cor_matrix = matrix()
result <- WGCNA::cor(t(module_agg), t(ann_qselect), method = "pearson")
pheatmap(result, color = colorRampPalette(c("navy", "white", "red"))(30), cluster_cols = FALSE)

# Calculate correlations with all filters
#Import the influence per stage
setwd("~/Oxford 2.0/Scripts/CNN_project/Data/motif_original_random")
table_target <- read.table("table_target.txt", header = FALSE)
setwd("~/Oxford 2.0/Scripts/CNN_project/Data/homer_comparison/")
table_target <- table_target[,-2]
table_target <- tidyr::spread(table_target, V3, V4)
colnames(table_target)[1] <- "Query.ID"
table_target$Query.ID <- paste0("filter",table_target$Query.ID)
table_target <- table_target[, c(1,6,3,9,8,7,5,4,2)]
rownames(table_target) <- table_target$Query.ID
table_target$Query.ID <- NULL

cor_matrix = matrix()
result_all <- WGCNA::cor(t(module_agg), t(table_target), method = "pearson")
result_all %>% select_if(~sum(!is.na(.)) > 0)

pheatmap(result_all, color = colorRampPalette(c("navy", "white", "red"))(30), cluster_cols = FALSE)


#Redo, with filter names only (so remove TF names), this means you need a smaller table first:
ann_slimmed <- sapply(strsplit(row.names(ann_qselect), ' '), function(x) x[1])
ann_qselect$row <- ann_slimmed
ann_slimdf <- ann_qselect[!duplicated(ann_slimmed),]
row.names(ann_slimdf) <- sapply(strsplit(row.names(ann_slimdf), ' '), function(x) x[1])
ann_slimdf$row <- NULL
ann_qselect$row <- NULL

#Now, calculate the correlations and plot:
cor_matrix_slim = matrix()
result_slim <- WGCNA::cor(t(module_agg), t(ann_slimdf), method = "pearson")
pheatmap(result_slim, color = colorRampPalette(c("navy", "white", "red"))(30), cluster_cols = FALSE)

#Redo, remove irrelevant modules
module_agg_select <- module_agg[c("black","green","red","turquoise"),]
cor_matrix_select = matrix()
result_select <- WGCNA::cor(t(module_agg_select), t(ann_slimdf), method = "pearson")
all_select <- WGCNA::cor(t(module_agg_select), t(table_target), method = "pearson")

pheatmap(result_select, color = colorRampPalette(c("navy", "white", "red"))(30), cluster_cols = FALSE, cluster_rows = FALSE)
pheatmap(all_select, color = colorRampPalette(c("navy", "white", "red"))(30), cluster_cols = FALSE, cluster_rows = FALSE)

### To make lineplots of overlap, you need identical stage names, for the CNNs, convert stage ###

### names to 2-letter names ###
colnames(ann_qselect)[3] <- "GT"
colnames(ann_qselect)[4] <- "PF"
new_lp <- rbind(ann_qselect[row.names(ann_qselect) == "filter26 RFX5_1",],
                module_agg[row.names(module_agg) == "turquoise",])
new_lp <- as.data.frame(t(new_lp))
new_lp$stage<- rownames(new_lp)
new_lp <- melt(new_lp)
#This is how to rescale, adjust accordingly. 25 is somewhat arbitrary
new_lp$value[new_lp$variable == "filter26 RFX5_1"] <- new_lp$value[new_lp$variable == 
                                                                     "filter26 RFX5_1"]*25
### Suggestion. Perhaps it is better in the long term to adjust accordingly by scaling to 1. ###
### Consider writing a script like that ###

#Once adjusted, make lineplot
p<-ggplot(new_lp, aes(x=stage, y=value, group=variable)) +
  geom_line(aes(color=variable), size =1)+
  geom_point(aes(color=variable), size =2.5) +
  ylim(-1.5,1.5)
p <- p + scale_color_brewer(palette="Set2")+
  labs(title="Filter influence per stage",x="Stage", y = "Influence") +
  theme_minimal(base_size = 22)
p + theme(legend.position="top") + theme(plot.title = element_text(hjust = 0.5))

# Make a plot with just TF family name (clip off the numbers)

# To quantify overlap: Merge the module table with the HOMER TF table. 
# Start with smaller motif_table.
# Keep essential columns
homer_foverlap <- as.data.frame(cbind(as.character(motif_table$module), as.character(motif_table$best_match)))
ann_foverlap1 <- sapply(strsplit(row.names(ann_qselect), ' '), function(x) x[1])
ann_foverlap2 <- sapply(strsplit(row.names(ann_qselect), ' '), function(x) x[2])
ann_fulloverlap <- as.data.frame(cbind(ann_foverlap1,ann_foverlap2))
homer_foverlap$V2 <- as.character(homer_foverlap$V2)
homer_foverlap$V2 <- sapply(strsplit(homer_foverlap$V2, '.motif'), function(x) x[1])
colnames(ann_fulloverlap) <- c("filter", "motif")
colnames(homer_foverlap) <- c("module", "motif")
mod_ann_overlap <- merge(ann_fulloverlap, homer_foverlap, by.x = "motif", by.y = "motif", all.x = TRUE, all.y = TRUE, 
      no.dups = TRUE)

# Make an empty matrix based on unique module and filter names
agg_mod_count <- aggregate(x = mod_ann_overlap$module, 
                           by = list(unique.values = mod_ann_overlap$module), 
                           FUN = length)
overlap_mat <- matrix(0, nrow = dim(homer_foverlap[!duplicated(homer_foverlap$module), ])[1], 
       ncol = dim(ann_fulloverlap[!duplicated(ann_fulloverlap$filter), ])[1])
rownames(overlap_mat) <- unique(homer_foverlap$module)
colnames(overlap_mat) <- unique(ann_fulloverlap$filter)

for (m in rownames(overlap_mat)) {
  for (n in colnames(overlap_mat)) {
  z <- mod_ann_overlap$module == m & mod_ann_overlap$filter == n
  overlap_mat[m,n] <- sum(z, na.rm=TRUE)
  }
}

#Redo for larger table
lhomer_foverlap <- as.data.frame(cbind(as.character(motif_matches$module), as.character(motif_matches$ID)))
lhomer_foverlap <- lhomer_foverlap[-1, ]
ann_foverlap1 <- sapply(strsplit(row.names(ann_qselect), ' '), function(x) x[1])
ann_foverlap2 <- sapply(strsplit(row.names(ann_qselect), ' '), function(x) x[2])
ann_fulloverlap <- as.data.frame(cbind(ann_foverlap1,ann_foverlap2))
lhomer_foverlap$V2 <- as.character(lhomer_foverlap$V2)
lhomer_foverlap$V2 <- sapply(strsplit(lhomer_foverlap$V2, '.motif'), function(x) x[1])
colnames(ann_fulloverlap) <- c("filter", "motif")
colnames(lhomer_foverlap) <- c("module", "motif")
lmod_ann_overlap <- merge(ann_fulloverlap, lhomer_foverlap, by.x = "motif", by.y = "motif", all.x = TRUE, all.y = TRUE, 
                         no.dups = TRUE)
write.table(lmod_ann_overlap, file = "mod_filter_overlap.txt")

lagg_mod_count <- aggregate(x = lmod_ann_overlap$module, 
                           by = list(unique.values = lmod_ann_overlap$module), 
                           FUN = length)

loverlap_mat <- matrix(0, nrow = dim(lhomer_foverlap[!duplicated(lhomer_foverlap$module), ])[1], 
                      ncol = dim(ann_fulloverlap[!duplicated(ann_fulloverlap$filter), ])[1])
rownames(loverlap_mat) <- unique(lhomer_foverlap$module)
colnames(loverlap_mat) <- unique(ann_fulloverlap$filter)

for (m in rownames(loverlap_mat)) {
  for (n in colnames(loverlap_mat)) {
    z <- lmod_ann_overlap$module == m & lmod_ann_overlap$filter == n
    loverlap_mat[m,n] <- sum(z, na.rm=TRUE)
  }
}

lann_homer_overlap <- t(apply(loverlap_mat, 1, "/", lagg_mod_count$x))
pheatmap(lann_homer_overlap, color = colorRampPalette(c("navy", "white", "red"))(30), cluster_cols = FALSE, cluster_rows = FALSE)

