#!/bin/bash
#SBATCH --job-name=phyml_windows_pseudoHap_n82_fossil_withOutgroup     # Job name
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --nodes=1                      # Run all processes on a single node
#SBATCH --ntasks=10                     # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G               # Memory required per CPU
#SBATCH --time=24:00:00                # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/phyml/BWD_HiC/n82_fossil/withOutgroup_all/%x_%j.out    # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/phyml/BWD_HiC/n82_fossil/withOutgroup_all/%x_%j.err     # Standard error log

# The %j in the --output & --error lines tells SLURM to substitute the job ID in the name of the output file
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules and set up environment variables
#module load phyml/3.3.20211231 # Not yet reinstalled on Mjolnir, so I installed it (v3.3.20220408) in my home dir and added it to my path.
## Source .bashrc to make sure phyml is in the PATH
source /home/KUID/.bashrc
module load python/2.7.18
MARTIN=/home/KUID/programs/genomics_general_simon_martin/phylo # Note: genomics.py must be in this directory, so copy it over before starting
DATA=/projects/lorenzen/people/KUID/wildebeest/admixture_graph/BWD_HiC # Where the haplo files are
OUTDIR=/projects/lorenzen/people/KUID/wildebeest/phyml/BWD_HiC/n82_fossil/withOutgroup_all
ASTRAL=/home/KUID/programs/Astral

# Navigate to working directory
cd $OUTDIR

# Use Simon Martin's phyml windows script to make ML trees in windows. Note: requires python2 and phyml in your PATH
$MARTIN/phyml_sliding_windows.py --phyml phyml \
 --genoFile $DATA/n82_fossil_BWD_HiC_autosomes.haplo \
 --prefix n82_fossil_BWD_HiC_autosomes_phyml_wOutgroup \
 --windSize 50000 --stepSize 1000000 \
 --minSites 20 --minPerInd 10 --optimise tlr \
 --threads 10 \
 --indFile n82_fossil_all.txt

# Use ASTRAL to get species tree and support for branches

## Filter out windows for which there are no trees (NA)
zcat n82_fossil_BWD_HiC_autosomes_phyml_wOutgroup.trees.gz | grep -v NA > n82_fossil.wOutgroup.50kb.trees

## Run ASTRAL with mapping file (population monophyly enforced)
java -jar $ASTRAL/astral.5.7.8.jar -i n82_fossil.wOutgroup.50kb.trees -a n82_fossil_mapping_all.txt -o n82_fossil.wOutgroup.50kb.ASTRAL --branch-annotate 32 2> n82_fossil.wOutgroup.50kb.ASTRAL.log
