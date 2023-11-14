#!/bin/bash
module load anaconda/4.10.3 
source activate erisONE_python_3.10_R_4.2.0
# module load Anaconda3/5.2.0
# source activate R_4.1.0_python_3.8.2


# load in the parameters
param_line=`cat missing_jobs.txt | head -n $LSB_JOBINDEX | tail -n1`
PARAM1=`echo $param_line | awk '{ print $1 }'`
PARAM2=`echo $param_line | awk '{ print $2 }'`
PARAM3=`echo $param_line | awk '{ print $3 }'`
PARAM4=`echo $param_line | awk '{ print $4 }'`
PARAM5=`echo $param_line | awk '{ print $5 }'`
PARAM6=`echo $param_line | awk '{ print $6 }'`

Rscript downsample.r $PARAM1 $PARAM2 $PARAM3 $PARAM4 $PARAM5 $PARAM6
echo $PARAM1 $PARAM2 $PARAM3 $PARAM4 $PARAM5 $PARAM6

# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=40000]' -W 120:00 -J cpca[51-75] "bash FINAL_downsample.sh"

# bsub -o OUT.%J.%I.out -q rerunnable -R 'rusage[mem=10000]' -W 1:00 -J umap[1-25] "bash FINAL_downsample.sh" 

# bsub -o OUT.%J.%I.out -q rerunnable -R 'rusage[mem=12000]' -W 48:00 -J pca[26-50] "bash FINAL_downsample.sh"

# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=100000]' -W 192:00 -J nmf[76-83] "bash FINAL_downsample.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=55000]' -W 120:00 -J nmf[84-93] "bash FINAL_downsample.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=55000]' -W 120:00 -J nmf[109-118] "bash FINAL_downsample.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=40000]' -W 48:00 -J nmf[94-108] "bash FINAL_downsample.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=40000]' -W 48:00 -J nmf[119-133] "bash FINAL_downsample.sh"








# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=120000]' -W 192:00 -J nmf[393-436] "bash FINAL_downsample.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=100000]' -W 120:00 -J nmf[467-496] "bash FINAL_downsample.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=100000]' -W 120:00 -J nmf[437-466] "bash FINAL_downsample.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=37000]' -W 48:00 -J pca[323-352] "bash FINAL_downsample.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=37000]' -W 48:00 -J umap[293-322] "bash FINAL_downsample.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=55000]' -W 120:00 -J cpca[353-392] "bash FINAL_downsample.sh"



# bsub -o OUT.%J.%I.out -q bigmem -R 'rusage[mem=90000]' -W 192:00 -J pca[232-301] "bash FINAL_downsample.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=90000]' -W 150:00 -J pca[232] "bash FINAL_downsample.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=90000]' -W 150:00 -J pca[262] "bash FINAL_downsample.sh"

# bsub -o OUT.%J.%I.out -q rerunnable -R 'rusage[mem=10000]' -W 1:00 -J umap[350-458] "bash FINAL_downsample.sh" 
# bsub -o OUT.%J.%I.out -q rerunnable -R 'rusage[mem=12000]' -W 12:00 -J cpca[1-108] "bash FINAL_downsample.sh" 




# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=30000]' -W 12:00 -J c_pca[343-512] "bash FINAL_downsample.sh" 

# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=35000]' -W 30:00 -J nmf[783-817] "bash FINAL_downsample.sh" 
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=35000]' -W 45:00 -J nmf[538-562] "bash FINAL_downsample.sh" 
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=35000]' -W 45:00 -J nmf[638-662] "bash FINAL_downsample.sh"
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=35000]' -W 30:00 -J nmf[728-747] "bash FINAL_downsample.sh"
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=40000]' -W 48:00 -J nmf[713-727] "bash FINAL_downsample.sh"
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=40000]' -W 45:00 -J nmf[513-537] "bash FINAL_downsample.sh" 
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=40000]' -W 45:00 -J nmf[613-637] "bash FINAL_downsample.sh" 
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=50000]' -W 48:00 -J nmf[748-762] "bash FINAL_downsample.sh"


 
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=45000]' -W 48:00 -J nmf[588-612] "bash FINAL_downsample.sh"
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=45000]' -W 48:00 -J nmf[688-712] "bash FINAL_downsample.sh"
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=50000]' -W 56:00 -J nmf[818-852] "bash FINAL_downsample.sh" 
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=85000]' -W 192:00 -J nmf[563-587] "bash FINAL_downsample.sh"
# bsub -o OUT.%J.%I.out -q normal -R 'rusage[mem=85000]' -W 192:00 -J nmf[663-687] "bash FINAL_downsample.sh"

 



