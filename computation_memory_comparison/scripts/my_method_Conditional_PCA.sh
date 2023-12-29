#conda activate CondPCA_eris2
module load anaconda/4.8.2
source activate scanpy_eris1 

# load in the parameters
param_line=`cat make_datasets.par | head -n $LSB_JOBINDEX | tail -n1`
CELLS=`echo $param_line | awk '{ print $1 }'`
GENES=`echo $param_line | awk '{ print $2 }'`

# Run the command and capture output and timing information
job_stats=$( { /usr/bin/time -v python my_method_Conditional_PCA.py $CELLS $GENES; } 2>&1 )
echo $job_stats
# extract items related to time and format (elapsed time isn't in seconds)
user_time_seconds=$(echo "$job_stats" | grep "User time (seconds):" | awk '{print $4}')
system_time_seconds=$(echo "$job_stats" | grep "System time (seconds):" | awk '{print $4}')
elapsed_time_seconds=$(echo "$job_stats" | grep "Elapsed (wall clock) time (h:mm:ss or m:ss):" | awk '{split($8, a, ":"); if (split($8, b, ":") == 3) print a[1]*3600 + a[2]*60 + a[3]; else if (split($8, b, ":") == 2) print a[1]*60 + a[2]}')
raw_elapsed_time=$(echo "$job_stats" | grep "Elapsed (wall clock) time (h:mm:ss or m:ss):" |awk '{print $8}')

# extract items related to memory
max_resident_set_size=$(echo "$job_stats" | grep "Maximum resident set size (kbytes):" | awk '{print $6}')
avg_resident_set_size=$(echo "$job_stats" | grep "Average resident set size (kbytes):" | awk '{print $6}')

#retrieve memory from r or python
memory_used_bytes=$(echo "$job_stats" | grep -o 'memory_used_bytes: *[0-9]*' | cut -d' ' -f2)

# add everything to output file
echo $CELLS $GENES Conditional my_method $user_time_seconds $system_time_seconds $elapsed_time_seconds $raw_elapsed_time $max_resident_set_size $avg_resident_set_size $memory_used_bytes >> output.txt
# bsub -o OUT.%J.%I.out -q rerunnable -R 'rusage[mem=12000]' -W 12:00 -J cmy[1-15] "bash my_method_Conditional_PCA.sh"
