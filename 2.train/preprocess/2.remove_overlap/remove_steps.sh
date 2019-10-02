bedtools intersect -a ER.bed -b original_islets.bed -v > negative_original_noverlap.bed
bedtools intersect -a ER.bed -b 1CPM_islets.bed -v > negative_1CPM_noverlap.bed
cp negative_1CPM_noverlap.bed ../final_steps/
cp negative_original_noverlap.bed ../final_steps/
rm negative_*
