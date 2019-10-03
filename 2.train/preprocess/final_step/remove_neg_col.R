#!/apps/well/R/3.4.3/bin/Rscript

test <- read.table("1CPM_islets_act.txt")
test <- test[,-9]
test$Neg <- rep(0, nrow(test))
write.table(test, file = "8colneg1CPM_islets_act.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
