# Use the shuffle function 
shuf -n 100000 negative_original_noverlap.bed > neg_original_noverlap_shuf.bed
shuf -n 100000 negative_1CPM_noverlap.bed > neg_1CPM_noverlap_shuf.bed
./column2zeros.R
