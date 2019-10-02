cd /well/mccarthy/users/maxlouis/oxford2/CNN_project/preprocessing/negative_set/data/2_remove_overlap
bedtools intersect -a ER.bed -b original_islets.bed -v > negative_original_noverlap.bed
bedtools intersect -a ER.bed -b 1CPM_islets.bed -v > negative_1CPM_noverlap.bed
cp negative_1CPM_noverlap.bed ../3_shuffle_and_zeros/
cp negative_original_noverlap.bed ../3_shuffle_and_zeros/
rm negative_*
