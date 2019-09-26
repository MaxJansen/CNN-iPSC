# the listin table target. 

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


