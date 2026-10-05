#!/bin/bash
#SBATCH --job-name=mapD_SRR27955277_trim_Mut_BWD_HiC     # Job name
#SBATCH --nodes=1                      # Run all processes on a single node
#SBATCH --ntasks=2                     # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G               # Memory required per CPU
#SBATCH --time=24:00:00                # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/paleomix/competitive_mapping/SRR27955277/BWD_HiC/mapDamage/mapD_SRR27955277_trim_Mut_%j.out    # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/paleomix/competitive_mapping/SRR27955277/BWD_HiC/mapDamage/mapD_SRR27955277_trim_Mut_%j.err     # Standard error log

# The %j in the --output & --error lines tells SLURM to substitute the job ID in the name of the output file
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used 
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules
module purge
module load mapdamage2/2.2.1

# Navigate to working directory and set up variables

#--- BWD_HiC-Human bam:
BAM=/projects/lorenzen/people/KUID/wildebeest/paleomix/competitive_mapping/SRR27955277/BWD_HiC
OUTDIR=/projects/lorenzen/people/KUID/wildebeest/paleomix/competitive_mapping/SRR27955277/BWD_HiC/mapDamage
REF="/projects/lorenzen/people/KUID/wildebeest/refgenomes/BWD_HiC.chr_Human_fossilMito.fasta"
### Insert commands
cd $OUTDIR
mapDamage -i $BAM/SRR27955277_trim_Mut.BWD_HiC-Human.bam -r $REF -d $OUTDIR --verbose --merge-reference-sequences --rescale --single-stranded -Q 20

