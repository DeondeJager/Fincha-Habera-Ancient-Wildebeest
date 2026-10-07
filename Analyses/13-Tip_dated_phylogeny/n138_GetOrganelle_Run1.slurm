#!/bin/bash
#SBATCH --job-name=beast_R1          # Job name
#SBATCH --nodes=1                       # Run all processes on a single node
#SBATCH --ntasks=1                      # Run on a single task
#SBATCH --cpus-per-task=4               # Number of CPU cores per task
#SBATCH --mem-per-cpu=20G               # Memory required per CPU
#SBATCH --time=48:00:00                  # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/beast/GetOrganelle/n138_tipDate/Run1/Run1_screenOut_%j.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/beast/GetOrganelle/n138_tipDate/Run1/Run1_screenErr_%j.err         # Standard error log

# The %j in the --output & --error lines tells SLURM to substitute the job ID in the name of the output file
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts 
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Activate conda environment/load modules
BEAST=/home/KUID/programs/beast/v2.6.7/bin

# Navigate to working directory
cd /projects/lorenzen/people/KUID/wildebeest/beast/GetOrganelle/n138_tipDate/Run1

# Insert commands here:
$BEAST/beast -threads 4 -overwrite n138_R1.xml

