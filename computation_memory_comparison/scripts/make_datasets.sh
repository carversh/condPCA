#!/bin/bash
module load anaconda/4.10.3
source activate sim
#module load Anaconda3/5.2.0
#source activate R_4.1.0_python_3.8.2


# load in the parameters
param_line=`cat make_datasets.par | head -n $LSB_JOBINDEX | tail -n1`
CELLS=`echo $param_line | awk '{ print $1 }'`
GENES=`echo $param_line | awk '{ print $2 }'`

Rscript make_datasets.R $CELLS $GENES

# bsub -o OUT.%J.%I.out -q rerunnable -R 'rusage[mem=12000]' -W 24:00 -J pipe[1-20] "bash make_datasets.sh"
