##################################################################
### This submits scripts in Basset_train to the cluster        ###
##################################################################

#!/bin/bash

#1CPM random
qsub -q gpu8.q -l gpu=2 -l gputype=p100  -V -cwd -N one_CPM_random_train -e 1CPM_random_train.err -o 1CPM_random_train.out Basset_train.sh
