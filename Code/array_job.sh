#!/bin/bash 
#SBATCH --job-name=inter-div # Job name 
#SBATCH --mail-type=BEGIN,END,FAIL 	     # Mail events (NONE, BEGIN, END, FAIL, ALL) 
#SBATCH --mail-user=wentao.yu@idiv.de # Where to send mail 
#SBATCH --chdir=/work/yuw/inter-div    # Set the working directory 
#SBATCH --time=7-12:00:00 		# Time limit day-hrs:min:sec 
#SBATCH --mem-per-cpu=15G  
#SBATCH --output=/work/%u/%x-%A-%a.out 	 # Standard output 
#SBATCH --error=/work/%u/%x-%A-%a.err    # error log  
	
module load GCC/12.2.0 OpenMPI/4.1.4 

module load R/4.2.2 # load the latest version of R

export MC_CORES=${SLURM_CPUS_PER_TASK:-9}

params="$1"
array_or_job_id=${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}

Rscript --vanilla  /work/yuw/inter-div/model_array.R "$params"
