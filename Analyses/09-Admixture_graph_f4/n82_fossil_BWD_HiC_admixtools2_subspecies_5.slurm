#!/bin/bash
#SBATCH --job-name=admixtools2_n82_fossil_BWD_HiC_subspecies_5     # Job name
#SBATCH --exclude=mjolnircomp05fl,mjolnircomp06fl
#SBATCH --nodes=1                      # Run all processes on a single node
#SBATCH --ntasks=20                     # Number of CPU cores per task
#SBATCH --mem-per-cpu=8G               # Memory required per CPU
#SBATCH --time=14-00:00:00                # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/admixture_graph/BWD_HiC/admixtools2/subspecies/%x_%j.out    # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/admixture_graph/BWD_HiC/admixtools2/subspecies/%x_%j.err     # Standard error log

# The %j in the --output & --error lines tells SLURM to substitute the job ID in the name of the output file
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used 
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

## Use locally installed conda environment of admixtools2
### Initiate conda environment on the compute node first
source /home/KUID/.bashrc
### Activate package
conda deactivate
conda activate admixtools2

# Navigate to working directory and set up variables
cd /projects/lorenzen/people/KUID/wildebeest/admixture_graph/BWD_HiC/admixtools2/subspecies

# Run R script
Rscript admixtools2_qpGraph_subspecies_5.R
rm Rplots.pdf

