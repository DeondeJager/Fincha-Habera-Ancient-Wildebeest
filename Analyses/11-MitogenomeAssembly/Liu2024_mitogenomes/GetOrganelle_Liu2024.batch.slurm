#!/bin/bash
#SBATCH --job-name=GetOrganelle_Liu2024         # Job name
#SBATCH --exclude=mjolnircomp05fl,mjolnircomp06fl # Exclude older nodes that might give memory problems (512GB)
#SBATCH --array=1-133%20		             # The range of indices for batch array (how many jobs)%No. of jobs to run in parallel.
#SBATCH --ntasks=1                      # Run on a single task
#SBATCH --cpus-per-task=4               # Number of CPU cores per task
#SBATCH --mem-per-cpu=100G                # Memory required per CPU
#SBATCH --time=24:00:00                 # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/GetOrganelle/GetOrganelle_Liu2024_%A_%a.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/GetOrganelle/GetOrganelle_Liu2024_%A_%a.err         # Standard error log

echo "Processing input file Liu_${SLURM_ARRAY_TASK_ID}.sh"
bash /home/KUID/scripts/GetOrganelle/Liu2024/Liu_${SLURM_ARRAY_TASK_ID}.sh
