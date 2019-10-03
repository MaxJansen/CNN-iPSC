#!/apps/well/R/3.4.3/bin/Rscript

cpm <- read.table("neg_1CPM_noverlap_shuf.bed", header = F)
original <- read.table("neg_original_noverlap_shuf.bed", header = F)

cpm$V7 <- rep(0, nrow(cpm))
original$V7 <- rep(0, nrow(original))

write.table(
  x = cpm,
  file = "neg_1CPM_zeros_shuf.bed",
  row.names = F,
  col.names = F,
  quote = F,
  sep = "\t"
)
write.table(
  x = original,
  file = "neg_original_zeros_shuf.bed",
  row.names = F,
  col.names = F,
  quote = F,
  sep = "\t"
)
