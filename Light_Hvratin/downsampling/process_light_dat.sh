#!/bin/bash
module load anaconda/4.10.3 
source activate erisONE_python_3.10_R_4.2.0
#module load Anaconda3/5.2.0
#source activate R_4.1.0_python_3.8.2


# load in the parameters
param_line=`cat downsampling.par | head -n $LSB_JOBINDEX | tail -n1`
PARAM1=`echo $param_line | awk '{ print $1 }'`
PARAM2=`echo $param_line | awk '{ print $2 }'`
PARAM3=`echo $param_line | awk '{ print $3 }'`
PARAM4=`echo $param_line | awk '{ print $4 }'`
PARAM5=`echo $param_line | awk '{ print $5 }'`

Rscript final_downsample_analysis.r $PARAM1 $PARAM2 $PARAM3 $PARAM4 $PARAM5

# bsub -o OUT.%J.%I.out -q rerunnable -R 'rusage[mem=5000]' -W 5:00 -J pipe[1-360] "bash process_light_dat.sh"



# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=32000]' -W 72:00 -J pipe[1086-1105] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[1016-1045] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=25000]' -W 48:00 -J pipe[1106-1145] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=32000]' -W 48:00 -J pipe[1046-1085] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 120:00 -J pipe[966-1015] "bash process_light_dat.sh"


# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=32000]' -W 48:00 -J pipe[1596-1856] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[1183] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[1197] "bash process_light_dat.sh"

# bsub -o OUT.%J.%I.out -q rerunnable -R 'rusage[mem=12000]' -W 48:00 -J pipe[1136] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[966-1025] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=32000]' -W 48:00 -J pipe[1026-1055] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=32000]' -W 48:00 -J pipe[1086-1115] "bash process_light_dat.sh"


# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=30000]' -W 48:00 -J pipe[1056-1085] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q rerunnable -R 'rusage[mem=10000]' -W 48:00 -J pipe[1116-1145] "bash process_light_dat.sh"




# bsub -o OUT.%J.%I.out -q rerunnable -R 'rusage[mem=5000]' -W 0:30 -J pipe[784-963] "bash process_light_dat.sh"


# bsub -o OUT.%J.%I.out -q rerunnable -R 'rusage[mem=5000]' -W 1:00 -J pipe[1] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=15000]' -W 24:00 -J pipe[362] "bash process_light_dat.sh"

# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[452-481] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[572-601] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[633-662] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[693-712] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[753-772] "bash process_light_dat.sh"

# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[405] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[408] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[418] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[523] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[524] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[525] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[299] "bash process_light_dat.sh"

# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=40000]' -W 60:00 -J pipe[494-521] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=40000]' -W 60:00 -J pipe[493] "bash process_light_dat.sh"





# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[432-451] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[603-612] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[663-672] "bash process_light_dat.sh"
# bsub -o OUT.%J.%I.out -q big -R 'rusage[mem=38000]' -W 72:00 -J pipe[723-732] "bash process_light_dat.sh"



