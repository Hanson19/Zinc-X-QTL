#!/bin/bash
#SBATCH --job-name=zinc_qc_SNPs      # Name of job
#SBATCH --mail-type=BEGIN,END,FAIL   # Set when I get emails (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=khanson15@ku.edu     # Email address
#SBATCH --partition=sjmac            # Either sixhour or sjmac
#SBATCH --nodes=1                    # Number of nodes to run on
#SBATCH --ntasks=1                   # Number of tasks to run 
#SBATCH --cpus-per-task=4            # Number of CPUs per task (>1 = multithreaded)
#SBATCH --mem-per-cpu=2gb            # Memory request
#SBATCH --time=6:00:00               # Time limit in hrs:min:sec 
#SBATCH --output=./run_log_files/qc_SNPs_%j.log    # Standard output and error logs

module load R/3.6
dir1="process"
echo -ne "CHROM\tPOS\tREF\tALT" > $dir1/SNPs.txt
cat bams.txt | cut -f4 -d"/" | cut -f1 -d"." | awk '{printf("\t%s_R\t%s_A",$1,$1)}' >> $dir1/SNPs.txt
echo -ne "\tblah\n" >> $dir1/SNPs.txt

# add the frequencies by chromosome to the header file
cat $dir1/temp.chrX.txt >> $dir1/SNPs.txt
cat $dir1/temp.chr2L.txt >> $dir1/SNPs.txt
cat $dir1/temp.chr2R.txt >> $dir1/SNPs.txt
cat $dir1/temp.chr3L.txt >> $dir1/SNPs.txt
cat $dir1/temp.chr3R.txt >> $dir1/SNPs.txt
mv $dir1/SNPs.txt .

# Rscript scripts/QC.freqs.R input NumberOfFounders
Rscript scripts/QC.freqs.R "SNPs.txt" 17
