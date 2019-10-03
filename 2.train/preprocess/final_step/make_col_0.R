#!/apps/well/R/3.4.3/bin/Rscript
# This script adds a column of zeros where the neg column used to be

test <- read.table("1CPM_islets_act.txt")
test$Neg <- rep(0, nrow(test))
write.table(test, file = "n1CPM_islets_act_test.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

