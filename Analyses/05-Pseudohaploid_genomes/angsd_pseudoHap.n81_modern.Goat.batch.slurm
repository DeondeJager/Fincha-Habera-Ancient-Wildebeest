#!/bin/bash
#SBATCH --job-name=angsd_pseudoHap_Goat_n81_modern    # Job name
#SBATCH --exclude=mjolnircomp05fl,mjolnircomp06fl # Exclude these two nodes, as they gave memory errors last time (they're older and have less memory (512GB) than the other nodes (1024GB).
#SBATCH --array=1-29%10   # The range of indices for batch array (how many jobs)%No. of jobs to run in parallel.
#SBATCH --nodes=1                      # Run all processes on a single node
#SBATCH --ntasks=20                     # Number of CPU cores per task
#SBATCH --mem-per-cpu=20G               # Memory required per CPU
#SBATCH --time=4-00:00:00                # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/admixture_graph/Goat/angsd_pseudoHap_Goat_n81_modern_%A.%a.out    # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/admixture_graph/Goat/angsd_pseudoHap_Goat_n81_modern_%A.%a.err     # Standard error log

echo "Processing input file angsd_pseudoHap.n81_modern.Goat_${SLURM_ARRAY_TASK_ID}.sh"
bash /home/KUID/scripts/angsd/fossilWildebeest/Goat/pseudoHap/angsd_pseudoHap.n81_modern.Goat_${SLURM_ARRAY_TASK_ID}.sh

# In an interactive session, combine out files after all runs have finished as follows, keeping only the header of the first file:
## See: https://unix.stackexchange.com/questions/60577/concatenate-multiple-files-with-same-header
#for file in n81_modern_*.haplo.gz; do pigz -d $file; done
#awk 'FNR==1 && NR!=1 { while (/^chr/) getline; } 1 {print}' n81_modern_*.haplo > n81_modern_Goat_autosomes.haplo
#for file in n81_modern_NC*.haplo; do pigz -p 2 $file; done

