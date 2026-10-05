#!/bin/bash
#SBATCH --job-name=mapD_FossilWildebeest     # Job name
#SBATCH --nodes=1                      # Run all processes on a single node
#SBATCH --ntasks=2                     # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G               # Memory required per CPU
#SBATCH --time=24:00:00                # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/paleomix/competitive_mapping/BWD_HiC/mapDamage/mapD_FossilWildbeest_%j.out    # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/paleomix/competitive_mapping/BWD_HiC/mapDamage/mapD_FossilWildbeest_%j.err     # Standard error log

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

#--- BWD_HiC-Human bam:
BAM=/projects/lorenzen/people/KUID/wildebeest/paleomix/competitive_mapping/BWD_HiC
OUTDIR=/projects/lorenzen/people/KUID/wildebeest/paleomix/competitive_mapping/BWD_HiC/mapDamage
REF="/projects/lorenzen/people/KUID/wildebeest/refgenomes/BWD_HiC.chr_Human_fossilMito.fasta"
## Insert commands
cd $OUTDIR
mapDamage -i $BAM/Fossil936.BWD_HiC-Human.merged.bam -r $REF -d $OUTDIR --verbose --merge-reference-sequences --rescale --single-stranded -Q 20

#--- Goat bam:
BAM=/projects/lorenzen/people/KUID/wildebeest/paleomix/goat
OUTDIR=/projects/lorenzen/people/KUID/wildebeest/plaeomix/goat/mapDamage
REF="/projects/lorenzen/people/KUID/wildebeest/refgenomes/GCF_001704415.2_ARS1.2_genomic.fasta"
## Insert commands
cd $OUTDIR
mapDamage -i $BAM/Fossil936.Goat.merged.bam -r $REF -d $OUTDIR --verbose --merge-reference-sequences --rescale --single-stranded -Q 20
