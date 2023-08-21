#!/bin/bash
source activate R_4.1.0_python_3.8.2
#module load Anaconda3/5.2.0
#source activate R_4.1.0_python_3.8.2


# load in the parameters
param_line=`cat downsampling.par | head -n $LSB_JOBINDEX | tail -n1`
METHOD=`echo $param_line | awk '{ print $1 }'`
DATA_PARAM=`echo $param_line | awk '{ print $2 }'`
NUM_CELLS=`echo $param_line | awk '{ print $3 }'`
SEED=`echo $param_line | awk '{ print $4 }'`
DIM=`echo $param_line | awk '{ print $5 }'`
TOTAL_WORKERS=20
NUM_GENES=20000

Rscript final_downsample_analysis.r $METHOD $DATA_PARAM $NUM_CELLS $SEED $DIM

FILENAME=$METHOD\_$NUM_CELLS\_$SEED\_$DATA_PARAM\_$DIM\.tab
NAME=$METHOD\_cells_$NUM_CELLS\_seed_$SEED\_$DATA_PARAM\_rank_$DIM\ 
echo $NAME

conda deactivate
source activate cnmf_env

# NOT POINTING TO CORRECT PATH FILENAME JUST FILE

cnmf prepare --output-dir /data/gusev/SCRNA/HRVATIN/downsampling/cNMF/cNMF_output/QC --name $NAME -c /data/gusev/SCRNA/HRVATIN/$FILENAME -k $DIM --n-iter 200 --seed 14 --numgenes $NUM_GENES --total-workers $TOTAL_WORKERS 

echo $NAME >> cnmf_names.par

# create factor params
for (( c=1; c<=$TOTAL_WORKERS; c++ ))
do
    echo  $NAME $c >> cnmf_step2.par
done

rm /data/gusev/SCRNA/HRVATIN/$FILENAME


# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=1000]' -W 0:30 -J pipe[1860-1895] "bash cNMF_initial.sh"