#!/bin/bash
#SBATCH --job-name=zinc_bwa          # Name of job
#SBATCH --mail-type=BEGIN,END,FAIL   # Set when I get emails (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=khanson15@ku.edu # Email address change 
#SBATCH --partition=sjmac            # Either sixhour or sjmac
#SBATCH --nodes=1                    # Number of nodes to run on
#SBATCH --ntasks=1                   # Number of tasks to run (1 with >1 CPUS per task is multithreaded)
#SBATCH --cpus-per-task=4            # Number of CPUs per task (>1 = multithreaded)
#SBATCH --mem-per-cpu=4gb            # Memory request (default is 2Gb)
#SBATCH --time=12:00:00              # Time limit in hrs:min:sec (default for sjmac queue is 8 hrs)
#SBATCH --output=./run_log_files/zinc_bwa_%j.log    # Standard output and error logs
#SBATCH --array=1-29                 # Set number (N) of parallel jobs (1-N)

module load bwa/0.7.17
module load samtools/1.9
module load bcftools/1.9

ref="ref/dm6.fa"
dir1="process"
dir2="data/bam/pbam"
files="data/raw/readname.mapping.txt" #This is that 5 column list that made above.

shortname=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f1` 
F1=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f2`
F2=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f3`
R1=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f4`
R2=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f5`

bwa mem -t 4 -M $ref $F1 $R1 | samtools view -bS - > $dir1/$shortname.1.temp.bam
bwa mem -t 4 -M $ref $F2 $R2 | samtools view -bS - > $dir1/$shortname.2.temp.bam
samtools merge -f $dir1/$shortname.temp.bam $dir1/$shortname.1.temp.bam $dir1/$shortname.2.temp.bam 
samtools sort $dir1/$shortname.temp.bam -o $dir2/$shortname.bam
samtools index $dir2/$shortname.bam
