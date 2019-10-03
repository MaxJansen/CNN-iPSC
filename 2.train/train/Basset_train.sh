module load torch/20170221-p100-gcc5.4.0
module load python/2.7.10-gcc4.9.3
export BEDTOOLS=/apps/well/bedtools/2.24.0/
export PATH=$BEDTOOLS:$PATH
export BASSETDIR=/well/got2d/agata/Basset/
export PATH=$BASSETDIR/src:$PATH
export PYTHONPATH=$BASSETDIR/src:$PYTHONPATH
export LUA_PATH="$BASSETDIR/src/?.lua;$LUA_PATH"
export PATH=${PATH}:/well/got2d/agata/bin/weblogo/
export PATH=${PATH}:/apps/well/meme/4.11.2_2/bin

for i in {1..30}; do
	printf "CNN: $i\n"
 	basset_train.lua -cudnn -drop_rate -job ../filt21_params.txt -stagnant_t 10 -save original.iter$i /well/mccarthy/users/maxlouis/oxford2/CNN_project/preprocessing/negative_set/data/final_step/1CPM/real_random/1CPM_islets.h5 > log_original.iter$i 2>&1 </dev/null
done;
