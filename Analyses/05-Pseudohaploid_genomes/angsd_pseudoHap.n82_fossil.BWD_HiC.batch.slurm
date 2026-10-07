#!/bin/bash
#SBATCH --job-name=angsd_pseudoHap_BWD_HiC_n82_fossil    # Job name
#SBATCH --exclude=mjolnircomp05fl,mjolnircomp06fl # Exclude these two nodes, as they gave memory errors last time (they're older and have less memory (512GB) than the other nodes (1024GB).
#SBATCH --array=1,2,4-29%10   # The range of indices for batch array (how many jobs)%No. of jobs to run in parallel.
#SBATCH --nodes=1                      # Run all processes on a single node
#SBATCH --ntasks=20                     # Number of CPU cores per task
#SBATCH --mem-per-cpu=20G               # Memory required per CPU
#SBATCH --time=2-00:00:00                # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/admixture_graph/BWD_HiC/angsd_pseudoHap_BWD_HiC_n82_fossil_%A.%a.out    # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/admixture_graph/BWD_HiC/angsd_pseudoHap_BWD_HiC_n82_fossil_%A.%a.err     # Standard error log

echo "Processing input file angsd_pseudoHap.n82_fossil.BWD_HiC_${SLURM_ARRAY_TASK_ID}.sh"
bash /home/KUID/scripts/angsd/fossilWildebeest/BWD_HiC/pseudoHap/angsd_pseudoHap.n82_fossil.BWD_HiC_${SLURM_ARRAY_TASK_ID}.sh

# In an interactive session, combine out files after all runs have finished as follows, keeping only the header of the first file:
## See: https://unix.stackexchange.com/questions/60577/concatenate-multiple-files-with-same-header
#for file in n82_fossil_HiC_scaffold_*.haplo.gz; do gzip -d $file; done
## Since we don't have a haplo file for HiC_scaffold_3 (the X chr), we can't use a range of numbers containing 3 in the filename for the command below, so we have to make two ranges, then combine those two files afterwards:
#awk 'FNR==1 && NR!=1 { while (/^chr/) getline; } 1 {print}' n82_fossil_HiC_scaffold_{1..2}.haplo > n82_fossil_BWD_HiC_autosomes.haplo.temp1
#awk 'FNR==1 && NR!=1 { while (/^chr/) getline; } 1 {print}' n82_fossil_HiC_scaffold_{4..29}.haplo > n82_fossil_BWD_HiC_autosomes.haplo.temp2
#awk 'FNR==1 && NR!=1 { while (/^chr/) getline; } 1 {print}' n82_fossil_BWD_HiC_autosomes.haplo.temp{1..2} > n82_fossil_BWD_HiC_autosomes.haplo
#rm n82_fossil_BWD_HiC_autosomes.haplo.temp1 n82_fossil_BWD_HiC_autosomes.haplo.temp2
#for file in n82_HiC_scaffold_*.haplo; do gzip $file; done
#rm *haplo.gz
#gzip n82_fossil_BWD_HiC_autosomes.haplo
