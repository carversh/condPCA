#!/bin/bash
conda env list
module load Anaconda3/5.2.0
conda env list
source activate R_4.1.0_python_3.8.2

# load in the parameters
param_line=`cat /data/gusev/USERS/scarver/Light_Hvratin/downsampling.par | head -n $LSB_JOBINDEX | tail -n1`
METHOD=`echo $param_line | awk '{ print $1 }'`
DATA_PARAM=`echo $param_line | awk '{ print $2 }'`
NUM_CELLS=`echo $param_line | awk '{ print $3 }'`
SEED=`echo $param_line | awk '{ print $4 }'`
DIM=`echo $param_line | awk '{ print $5 }'`
NUM_GENES=`echo $param_line | awk '{ print $6 }'`
TOTAL_WORKERS=20



Rscript /data/gusev/USERS/scarver/Light_Hvratin/FINAL_downsample.r $METHOD $DATA_PARAM $NUM_CELLS $SEED $DIM $NUM_GENES

FILENAME=$METHOD\_$NUM_CELLS\_$SEED\_$DATA_PARAM\_$DIM\_$NUM_GENES\.tab
NAME=$METHOD\_cells_$NUM_CELLS\_seed_$SEED\_$DATA_PARAM\_rank_$DIM\_$NUM_GENES 
echo $NAME

conda deactivate
source activate cnmf_env


cnmf prepare --output-dir /data/gusev/USERS/scarver/Light_Hvratin/cNMF_final/cNMF_output/ --name $NAME -c /data/gusev/USERS/scarver/Light_Hvratin/cNMF_final/$FILENAME -k $DIM --n-iter 200 --seed 14 --numgenes $NUM_GENES --total-workers $TOTAL_WORKERS 



# create factor params
for (( c=1; c<=$TOTAL_WORKERS; c++ ))
do
    ((NUM=c-1))
    echo  $METHOD $DATA_PARAM $NUM_CELLS $SEED $DIM $NUM_GENES $NUM >> cnmf_step2.par
done

rm /data/gusev/USERS/scarver/Light_Hvratin/cNMF_final/$FILENAME

#echo $METHOD $DATA_PARAM $NUM_CELLS $SEED $DIM $NUM_GENES >> step1_finished.txt

# bsub -o OUT.%J.%I.out -q bigmem -R 'rusage[mem=33000]' -W 1:00 -J pipe[1902]  "bash cNMF_step1.sh"
