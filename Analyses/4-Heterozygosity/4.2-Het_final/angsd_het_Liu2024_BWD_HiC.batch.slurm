#!/bin/bash
#SBATCH --job-name=het_Liu2024_BWD_HiC			   # Job name
#SBATCH --exclude=mjolnircomp05fl,mjolnircomp06fl # Exclude these two nodes, as they gave memory errors last time (they're older and have less memory (512GB) than the other nodes (1024GB).
#SBATCH --array=1,2,3,4,5,6,8,9,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,33,34,35,36,37,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,74,89,90,99,100,101,102,103,104,105,106,107,108,111,121,126,127,128,131%10   # The range of indices for batch array (how many jobs)%No. of jobs to run in parallel.
#SBATCH --ntasks=1                         # Run on a single task
#SBATCH --cpus-per-task=10                 # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G                  # Memory required per CPU
#SBATCH --time=18:00:00                    # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/het/BWD_HiC/angsd_het_Liu2024_BWD_HiC_%A.%a.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/het/BWD_HiC/angsd_het_Liu2024_BWD_HiC_%A.%a.err         # Standard error log

echo "Processing input file angsd_het_Liu_${SLURM_ARRAY_TASK_ID}_BWD_HiC.sh"
bash /home/KUID/scripts/angsd/fossilWildebeest/BWD_HiC/het/angsd_het_Liu_${SLURM_ARRAY_TASK_ID}_BWD_HiC.sh

