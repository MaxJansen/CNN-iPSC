##### load torch modules & dependencies
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

# Run iter1 from the new 8 column negative set format interactively
nohup basset_test.lua -cudnn original.iter1_best.th /well/mccarthy/users/maxlouis/oxford2/CNN_project/preprocessing/negative_set/data/final_step/1CPM/real_random/1CPM_islets.h5 AUC_per_feature/iter1 &
