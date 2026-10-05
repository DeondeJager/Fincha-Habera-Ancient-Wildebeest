#!/bin/bash
#SBATCH --job-name=repeatmasker         # Job name
#SBATCH --nodes=1                       # Run all processes on a single node
#SBATCH --ntasks=1                      # Run on a single task
#SBATCH --cpus-per-task=16               # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G                # Memory required per CPU
#SBATCH --time=144:00:00                 # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/refgenomes/repeatmasker/repeatmasker_Goat_%j.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/refgenomes/repeatmasker/repeatmasker_Goat_%j.err         # Standard error log

# The %j in the --output & --error lines tells SLURM to substitute the job ID in the name of the output file
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts 
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules
module load repeatmasker/4.1.2.p1

# Set environment variablea and navigate to working directory
DATA=/projects/lorenzen/people/KUID/wildebeest/refgenomes
OUT=/projects/lorenzen/people/KUID/wildebeest/refgenomes/repeatmasker

# Run repeatmasker
cd $DATA
RepeatMasker -pa 4 \
 -dir $OUT \
 -species "Capra hircus" \
 $DATA/GCF_001704415.2_ARS1.2_genomic.fasta
