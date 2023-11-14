#!/bin/bash
module load Anaconda3/5.2.0
source activate cnmf_env

# load in the parameters
param_line=`cat step1.par | head -n $LSB_JOBINDEX | tail -n1`
METHOD=`echo $param_line | awk '{ print $1 }'`
DATA_PARAM=`echo $param_line | awk '{ print $2 }'`
NUM_CELLS=`echo $param_line | awk '{ print $3 }'`
SEED=`echo $param_line | awk '{ print $4 }'`
DIM=`echo $param_line | awk '{ print $5 }'`
NUM_GENES=`echo $param_line | awk '{ print $6 }'`

NAME=$METHOD\_cells_$NUM_CELLS\_seed_$SEED\_$DATA_PARAM\_rank_$DIM\_$NUM_GENES

cnmf combine --output-dir /PHShome/sv433/Gusev_Lab/Light_Hvratin/cNMF_final/cNMF_output/ --name $NAME
cnmf k_selection_plot --output-dir /PHShome/sv433/Gusev_Lab/Light_Hvratin/cNMF_final/cNMF_output/ --name $NAME

cnmf consensus --output-dir /PHShome/sv433/Gusev_Lab/Light_Hvratin/cNMF_final/cNMF_output/ --name $NAME --components $DIM --local-density-threshold 2.00 --show-clustering

conda deactivate 
source activate R_4.1.0_python_3.8.2

# read in the number of vargenes from the outputted cnmf file
VARGENES=$(head -1 /PHShome/sv433/Gusev_Lab/Light_Hvratin/cNMF_final/cNMF_output/$NAME\/$NAME\.spectra.k_$DIM\.dt_2_0.consensus.txt | awk '{ print NF}')

Rscript process_final_cnmf.r $NAME $METHOD $DATA_PARAM $NUM_CELLS $SEED $DIM $NUM_GENES

echo $METHOD $DATA_PARAM $NUM_CELLS $SEED $DIM $NUM_GENES >> completed_jobs_cnmf.txt

# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=25000]' -W 6:00 -J cnmf[1-34] "bash cNMF_step_final.sh"