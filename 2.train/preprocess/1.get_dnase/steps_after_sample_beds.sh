module load torch/20170221-k80-gcc5.4.0
module load python/2.7.10-gcc4.9.3
export BEDTOOLS=/apps/well/bedtools/2.24.0/
export PATH=$BEDTOOLS:$PATH
export BASSETDIR=/well/got2d/agata/Basset/
export PATH=$BASSETDIR/src:$PATH
export PYTHONPATH=$BASSETDIR/src:$PYTHONPATH
export LUA_PATH="$BASSETDIR/src/?.lua;$LUA_PATH"
export PATH=${PATH}:/well/got2d/agata/bin/weblogo/
export PATH=${PATH}:/apps/well/meme/4.11.2_2/bin

# Go to directory
cd /well/mccarthy/users/maxlouis/oxford2/CNN_project/preprocessing/negative_set/data/1_get_dnase

######################     create input files     #######################
### samples.txt - file with all the BED files to include in format: name \t filename
preprocess_features.py -y -m 200 -s 1000 -n -o ER -c /well/got2d/agata/Basset/data/genomes/human.hg19.genome ../1_get_dnase/sample_beds.txt
bedtools getfasta -fi /well/got2d/agata/Basset/data/genomes/hg19.fa -bed ER.bed -s -fo ER.fa
cp ER.bed ../2_remove_overlap/
