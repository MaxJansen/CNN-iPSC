# Load bedtools package:
export BEDTOOLS=/apps/well/bedtools/2.24.0/
export PATH=$BEDTOOLS:$PATH
# Go to the directory:
cd /well/mccarthy/users/maxlouis/oxford2/CNN_project/preprocessing/negative_set/data/2_remove_overlap

# Remove peaks occurring in samples from negative set
bedtools intersect -a ER.bed -b original_islets.bed -v > negative_original_noverlap.bed
bedtools intersect -a ER.bed -b 1CPM_islets.bed -v > negative_1CPM_noverlap.bed

# Move files for next step and remove clutter
cp negative_1CPM_noverlap.bed ../3_shuffle_and_zeros/
cp negative_original_noverlap.bed ../3_shuffle_and_zeros/
rm negative_*
