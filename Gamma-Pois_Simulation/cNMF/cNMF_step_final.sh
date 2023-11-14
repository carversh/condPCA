#!/bin/bash
module load Anaconda3/5.2.0
source activate cnmf_env

# load in the parameters
#param_line=`cat /PHShome/sv433/Gusev_Lab/Simulation/Gamma-Poisson/simulation.par | head -n $LSB_JOBINDEX | tail -n1`
param_line=`cat /PHShome/sv433/Gusev_Lab/Simulation/Gamma-Poisson/simulation.par | head -n $LSB_JOBINDEX | tail -n1`
SEED=`echo $param_line | awk '{ print $1 }'`
METHOD=`echo $param_line | awk '{ print $2 }'`
DATA_PARAM=`echo $param_line | awk '{ print $3 }'`
DIM=`echo $param_line | awk '{ print $4 }'`
NUM_CELLS=`echo $param_line | awk '{ print $5 }'`
PERC_GENES=`echo $param_line | awk '{ print $6 }'`
PERC_CELLS=`echo $param_line | awk '{ print $7 }'`
FLAG_100_CTS=`echo $param_line | awk '{ print $8 }'`
TOTAL_WORKERS=10
NUM_GENES=2527


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

cnmf combine --output-dir ./cNMF_output/ --name $NAME
cnmf k_selection_plot --output-dir ./cNMF_output/ --name $NAME

#cnmf consensus --output-dir ./cNMF_output/ --name $NAME --components $DIM --local-density-threshold 2.00 --show-clustering
cnmf consensus --output-dir ./cNMF_output/ --name $NAME --components $DIM --local-density-threshold 0.5 --show-clustering

conda deactivate 
source activate R_4.1.0_python_3.8.2

# read in the number of vargenes from the outputted cnmf file
#VARGENES=$(head -1 ./cNMF_output/$NAME\/$NAME\.spectra.k_$DIM\.dt_2_0.consensus.txt | awk '{ print NF}')
VARGENES=$(head -1 ./cNMF_output/$NAME\/$NAME\.spectra.k_$DIM\.dt_0_5.consensus.txt | awk '{ print NF}')

Rscript process_final_cnmf.r $SEED $METHOD $DATA_PARAM $DIM $NUM_CELLS $PERC_GENES $PERC_CELLS $NAME $FLAG_100_CTS 

echo $SEED $METHOD $DATA_PARAM $DIM $NUM_CELLS $PERC_GENES $PERC_CELLS $FLAG_100_CTS >> completed_jobs_cnmf.txt

# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=20000]' -W 6:00 -J fcnmf[641] "bash cNMF_step_final.sh"

# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=25000]' -W 2:00 -J fcnmf[1-70] "bash cNMF_step_final.sh"

# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=20000]' -W 7:00 -J fcnmf[120-390] "bash cNMF_step_final.sh"

#[1601-1920]

# bsub -o OUT.%J.%I.out -q bigmem -R 'rusage[mem=38000]' -W 12:00 -J s_cnmf[1601-1605] "bash cNMF_step_final.sh"
# bsub -o OUT.%J.%I.out -q bigmem -R 'rusage[mem=50000]' -W 8:00 -J s_cnmf[1606-1879] "bash cNMF_step_final.sh"
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=50000]' -W 8:00 -J s_cnmf[1880-1920] "bash cNMF_step_final.sh"