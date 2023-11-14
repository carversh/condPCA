#!/bin/bash
module load Anaconda3/5.2.0
source activate cnmf_env

# load in the parameters
param_line=`cat cnmf_step2.par | head -n $LSB_JOBINDEX | tail -n1`
SEED=`echo $param_line | awk '{ print $1 }'`
METHOD=`echo $param_line | awk '{ print $2 }'`
DATA_PARAM=`echo $param_line | awk '{ print $3 }'`
DIM=`echo $param_line | awk '{ print $4 }'`
NUM_CELLS=`echo $param_line | awk '{ print $5 }'`
PERC_GENES=`echo $param_line | awk '{ print $6 }'`
PERC_CELLS=`echo $param_line | awk '{ print $7 }'`
WORKER_INDEX=`echo $param_line | awk '{ print $8 }'`
FLAG_100_CTS=`echo $param_line | awk '{ print $9 }'`
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


cnmf factorize --output-dir ./cNMF_output/ --name $NAME --worker-index $WORKER_INDEX --total-workers 1

echo $SEED $METHOD $DATA_PARAM $DIM $NUM_CELLS $PERC_GENES $PERC_CELLS $WORKER_INDEX $FLAG_100_CTS >> step2_finished.txt

# bsub -o OUT.%J.%I.out -q bigmem -R 'rusage[mem=33000]' -W 9:30 -J scnmf2[1-50] "bash cNMF_step2.sh"
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=15000]' -W 9:30 -J scnmf2[51-320] "bash cNMF_step2.sh"
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=15000]' -W 24:00 -J scnmf2[321-626] "bash cNMF_step2.sh"
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=32000]' -W 6:30 -J scnmf2[322-418] "bash cNMF_step2.sh"
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=32000]' -W 6:30 -J scnmf2[419-460] "bash cNMF_step2.sh"

# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=32000]' -W 6:30 -J scnmf2[1-70] "bash cNMF_step2.sh"
#400 iter 
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=25000]' -W 24:00 -J scnmf2[2-280] "bash cNMF_step2.sh"



