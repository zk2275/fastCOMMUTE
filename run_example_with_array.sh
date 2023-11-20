#!/bin/bash
#SBATCH --account=free
#SBATCH -c 1                                                     # Number of cores (-c)
#SBATCH -t 1-00:00                                               # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --mem=10G #50G for original                                               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --job-name=fastCOMMUTE
#SBATCH -o real_data.out                   # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e real_data.err                   # File to which STDERR will be written, %j inserts jobid
#SBATCH --mail-user=zk2275@cumc.columbia.edu
#SBATCH --mail-type=END
#SBATCH --array=1-2 #1-200

#module load python/3.8.5-fasrc01 #Load Perl module
module load R/4.1.0-fasrc01 #Load R module
Rscript example_with_array_job.R 'test_item'