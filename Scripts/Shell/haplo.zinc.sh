
#!/bin/bash
#SBATCH --job-name=haplotyper        # Name of job
#SBATCH --mail-type=BEGIN,END,FAIL   # Set when I get emails (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=khanson15@ku.edu     # Email address
#SBATCH --partition=sjmac            # Either sixhour or sjmac
#SBATCH --nodes=1                    # Number of nodes to run on
#SBATCH --ntasks=1                   # Number of tasks to run 
#SBATCH --cpus-per-task=1            # Number of CPUs per task 
#SBATCH --mem-per-cpu=2gb            # Memory request 
#SBATCH --time=6:00:00               # Time limit in hrs:min:sec
#SBATCH --output=./run_log_files/haplotyper_%j.log    # Standard output and error logs
#SBATCH --array=1-29                 # Set number (N) of parallel jobs (1-N)

module load R/4.0
# TELL R where to find libraries
# because 'haplo.zinc.R' called below calls 'library(limSolve)'
export R_LIBS_USER=/panfs/pfs.local/work/sjmac/sjmac/R

file=$1
folder=$2
SNPtable=$3
foundernames=$4

pool=`head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1 | cut -f 1` 
OutName=`head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1 | cut -f 2`

echo "Rscript scripts/haplo.zinc.R $pool $OutName $folder $SNPtable $foundernames"
Rscript --verbose scripts/haplo.zinc.R $pool $OutName $folder $SNPtable $foundernames
