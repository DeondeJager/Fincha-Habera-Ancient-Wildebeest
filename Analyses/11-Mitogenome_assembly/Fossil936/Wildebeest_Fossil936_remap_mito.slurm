#!/bin/bash
#SBATCH --job-name=paleomix_Fossil936_remap_mito         # Job name
#SBATCH --nodes=1                       # Run all processes on a single node
#SBATCH --ntasks=1                      # Run on a single task
#SBATCH --cpus-per-task=2               # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G                # Memory required per CPU
#SBATCH --time=03:00:00                 # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/paleomix/936_final/%x_%j.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/paleomix/936_final/%x_%j.err         # Standard error log

# The %j in the --output & --error lines tells SLURM to substitute the job ID in the name of the output file
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts 
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script
hostname; date

echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules
module load paleomix/1.3.8

# Navigate to working directory where yaml file is
cd /projects/lorenzen/people/KUID/wildebeest/paleomix/936_final

# Execute paleomix
paleomix bam run 936_final_HiSeq4000.yaml --max-threads 2 --adapterremoval-max-threads 2 --bwa-max-threads 2

