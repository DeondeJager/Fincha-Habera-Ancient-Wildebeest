#!/bin/bash
#SBATCH --job-name=angsd_Abbababba_BWD_HiC_refBiasCheck     # Job name
#SBATCH --nodes=1                      # Run all processes on a single node
#SBATCH --ntasks=10                     # Number of CPU cores per task
#SBATCH --mem-per-cpu=20G               # Memory required per CPU
#SBATCH --time=04:00:00                # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/Dstats/angsd/Abbababba/BWD_HiC/%x_%j.out    # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/Dstats/angsd/Abbababba/BWD_HiC/%x_%j.err     # Standard error log

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
module load angsd/0.921
module load gcc/11.2.0 # For R
module load R/3.6.1 
RSCRIPTS=/home/KUID/programs/angsd/angsd-0.921/R

# Navigate to working directory and set up variables
OUTDIR=/projects/lorenzen/people/KUID/wildebeest/Dstats/angsd/Abbababba/BWD_HiC
REGIONS=/projects/lorenzen/people/KUID/wildebeest/refgenomes

#-- Run angsd abbababba test with damaged modern and undamaged counterpart, and one representative from each of the other subspecies (n=10)
## Note: min and max depths were calculated as the sum across all samples of the per-sample 0.33*median depth for min and median+2*SD depth for max.
## Note2: -setMaxDepth is for reading sites across the bam file, while -maxDepth is for the abbababba test to "avoid overflow of the counts". Not sure if these two will clash or not.
cd $OUTDIR
### One chromosome
angsd -doAbbababa 1 \
 -rmTrans 1 \
 -blockSize 5000000 \
 -bam bamlist_n10_refBiasCheck.txt \
 -r HiC_scaffold_1 \
 -doCounts 1 -setMinDepth 58 -setMaxDepth 349 \
 -baq 1 -ref $REGIONS/BWD_HiC.chr_Human_fossilMito.fasta \
 -maxDepth 349 \
 -out n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1 \
 -useLast 1 \
 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 \
 -nThreads 10

#-- Run jackknife script
Rscript $RSCRIPTS/jackKnife.R file=n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1.abbababa indNames=indNames_noOutgroup_refBiasCheck.txt boot=1 outfile=n10_BWD_HiC.refBiasCheck.rmTrans.noMaf.baq1.jackknife

