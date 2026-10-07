#!/bin/bash
#SBATCH --job-name=Fossil936_BWD_HiC_coverage_SeXY         # Job name
#SBATCH --nodes=1                       # Run all processes on a single node
#SBATCH --ntasks=1                      # Run on a single task
#SBATCH --cpus-per-task=2               # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G                # Memory required per CPU
#SBATCH --time=03:00:00                 # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/SeXY/BWD_HiC/%x_%j.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/SeXY/BWD_HiC/%x_%j.err         # Standard error log

# The %j in the --output & --error lines tells SLURM to substitute the job ID in the name of the output file
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script
hostname; date

echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Navigate to working directory
#cd /projects/lorenzen/people/KUID/wildebeest/SeXY

# Execute coverage script
# Calculate depths of bam files
## Randomly sample 1M sites 10x and calculate depth from said sites
## The Coverage_calculation.sh script requires samtools, which is loaded as a module inside that script

## $1 Autosomes.bed
## $2 X_trim.bed
## $3 Output_folder
## $4 Bamfile
## $5 Output_prefix

## For BWD_HiC version of Fossil936 (HiC_Scaffold_3 is X, everything else is autosome)
cd /projects/lorenzen/people/KUID/wildebeest/SeXY/BWD_HiC
bash ~/scripts/SeXY/CoverageCalculation.sh HiC_autosomes.bed HiC_scaffold_3.bed ./ Fossil936.BWD_HiC-Human.bam Fossil936.BWD_HiC-Human

# This produces the file "Fossil936.BWD_HiC-Human_ratios.txt" where the 10 chrX:Autosomes ratios can be found. This was plotted in R to produce Fig. 1C.

