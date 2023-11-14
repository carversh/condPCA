#!/bin/bash
module load Anaconda3/5.2.0
source activate R_4.1.0_python_3.8.2

# load in the parameters
# param_line=`cat /data/gusev/USERS/scarver/Gusev_Lab/Simulation/Gamma-Poisson/simulation.par | head -n $LSB_JOBINDEX | tail -n1`
param_line=`cat /data/gusev/USERS/scarver/Gusev_Lab/Simulation/Gamma-Poisson/missing_jobs.txt | head -n $LSB_JOBINDEX | tail -n1`
SEED=`echo $param_line | awk '{ print $1 }'`
METHOD=`echo $param_line | awk '{ print $2 }'`
DATA_PARAM=`echo $param_line | awk '{ print $3 }'`
DIM=`echo $param_line | awk '{ print $4 }'`
NUM_CELLS=`echo $param_line | awk '{ print $5 }'`
PERC_GENES=`echo $param_line | awk '{ print $6 }'`
PERC_CELLS=`echo $param_line | awk '{ print $7 }'`
FLAG_100_CTS=`echo $param_line | awk '{ print $8 }'`
TOTAL_WORKERS=1
NUM_GENES=2527


Rscript /data/gusev/USERS/scarver/Gusev_Lab/Simulation/Gamma-Poisson/Simulation.r $SEED $METHOD $DATA_PARAM $DIM $NUM_CELLS $PERC_GENES $PERC_CELLS $FLAG_100_CTS

# remove period becuase that messes up file structure
EDITED_perc_genes=$(echo $PERC_GENES | sed 's/\./p/g')
EDITED_perc_cells=$(echo $PERC_CELLS | sed 's/\./p/g')

# check if 100 ct flag is present
if [ ! -z "$FLAG_100_CTS" ] ; then
      FILENAME=$METHOD\_$NUM_CELLS\_$SEED\_$DATA_PARAM\_$DIM\_gene_$EDITED_perc_genes\_cell_$EDITED_perc_cells\_flag_$FLAG_100_CTS\.tab
      NAME=$METHOD\_cells_$NUM_CELLS\_seed_$SEED\_$DATA_PARAM\_rank_$DIM\_gene_$EDITED_perc_genes\_cell_$EDITED_perc_cells\_flag_$FLAG_100_CTS
      
else
    FILENAME=$METHOD\_$NUM_CELLS\_$SEED\_$DATA_PARAM\_$DIM\_gene_$EDITED_perc_genes\_cell_$EDITED_perc_cells\.tab
    NAME=$METHOD\_cells_$NUM_CELLS\_seed_$SEED\_$DATA_PARAM\_rank_$DIM\_gene_$EDITED_perc_genes\_cell_$EDITED_perc_cells
fi

echo $NAME

conda deactivate
source activate cnmf_env

# before, iterations was at 200 and done for 7 cells types
# iterations was changed to 400 for 100 cell types
cnmf prepare --output-dir ./cNMF_output/ --name $NAME -c ./$FILENAME -k $DIM --n-iter 400 --seed 14 --numgenes $NUM_GENES --total-workers $TOTAL_WORKERS 



# create factor params
for (( c=1; c<=$TOTAL_WORKERS; c++ ))
do
    ((NUM=c-1))
    echo  $SEED $METHOD $DATA_PARAM $DIM $NUM_CELLS $PERC_GENES $PERC_CELLS $NUM $FLAG_100_CTS >> cnmf_step2.par
done

rm ./$FILENAME

echo $SEED $METHOD $DATA_PARAM $DIM $NUM_CELLS $PERC_GENES $PERC_CELLS $FLAG_100_CTS >> step1_finished.txt

# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=30000]' -W 5:00 -J pipe[1-119]  "bash cNMF_step1.sh"
