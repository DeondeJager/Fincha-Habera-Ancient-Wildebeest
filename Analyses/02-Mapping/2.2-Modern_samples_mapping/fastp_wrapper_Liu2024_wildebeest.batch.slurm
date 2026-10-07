#!/bin/bash
#SBATCH --job-name=fastp_Liu2024_wildebeest         # Job name
#SBATCH --array=1-133%10                # The range of indices for batch array (how many jobs)%Maximum no. of jobs to run in parallel.
#SBATCH --ntasks=1                      # Run on a single task
#SBATCH --cpus-per-task=10               # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G                # Memory required per CPU
#SBATCH --time=03:00:00                 # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/clean_data/Liu_2024/fastp_%A.%a.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/clean_data/Liu_2024/fastp_%A.%a.err         # Standard error log

# In the --output & --error lines note the use of %A for the master job ID of the array, and the %a for the task ID in the output filename.
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script

echo "Processing input file Liu_${SLURM_ARRAY_TASK_ID}.sh"
bash /home/KUID/scripts/readsQC/Liu2024_fastp/Liu_${SLURM_ARRAY_TASK_ID}.sh

