#!/bin/bash
module load Anaconda3/5.2.0
source activate cnmf_env

# load in the parameters
param_line=`cat cnmf_step2.par | head -n $LSB_JOBINDEX | tail -n1`
METHOD=`echo $param_line | awk '{ print $1 }'`
DATA_PARAM=`echo $param_line | awk '{ print $2 }'`
NUM_CELLS=`echo $param_line | awk '{ print $3 }'`
SEED=`echo $param_line | awk '{ print $4 }'`
DIM=`echo $param_line | awk '{ print $5 }'`
NUM_GENES=`echo $param_line | awk '{ print $6 }'`
WORKER_INDEX=`echo $param_line | awk '{ print $7 }'`

NAME=$METHOD\_cells_$NUM_CELLS\_seed_$SEED\_$DATA_PARAM\_rank_$DIM\_$NUM_GENES
echo $NAME $WORKER_INDEX

cnmf factorize --output-dir /data/gusev/USERS/scarver/Light_Hvratin/cNMF_final/cNMF_output --name $NAME --worker-index $WORKER_INDEX --total-workers 20

echo $METHOD $DATA_PARAM $NUM_CELLS $SEED $DIM $NUM_GENES $WORKER_INDEX >> step2_finished.txt

# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=25000]' -W 1:00 -J 2lcnmf[1-100] "bash cNMF_step2.sh"
# bsub -o OUT.%J.%I.out -q bigmem -R 'rusage[mem=32000]' -W 24:00 -J pipe[1-2720] "bash cNMF_step2.sh"

# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=20000]' -W 24:00 -J lscnf2[1-500] "bash cNMF_step2.sh"


