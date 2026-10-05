#!/bin/bash
#SBATCH --job-name=Subset_Goat              # Job name
#SBATCH --nodes=1                       # Run on a single node
#SBATCH --ntasks=1                      # Run on a single task
#SBATCH --cpus-per-task=2              # Number of CPU cores per task
#SBATCH --mem-per-cpu=30G               # Memory required per CPU
#SBATCH --time=02:00:00                 # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/paleomix/goat/mapDamage/subset_%j.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/paleomix/goat/mapDamage/subset_%j.err         # Standard error log

# The %j in the --output & --error lines tells SLURM to substitute the job ID in the name of the output file
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

#--- Subset goat bam to autosomes and with repeat regions masked ---#

# Load modules
module purge
module load bedtools/2.31.0
module load samtools/1.15
module load gatk4/4.3.0.0

#--- Goat reference genome ---#
REGIONS="/projects/lorenzen/people/KUID/wildebeest/refgenomes/repeatmasker/GCF_001704415.2_ARS1.2_genomic_rm.sorted.complement.noX.bed"
GENOMEFILE="/projects/lorenzen/people/KUID/wildebeest/refgenomes/GCF_001704415.2_ARS1.2_genomic.chr.sizes.txt"
OUTDIR=/projects/lorenzen/people/KUID/wildebeest/paleomix/goat/mapDamage

cd $OUTDIR

# Subset bam file using regions file
date
echo "Subsetting bam file..."
bedtools intersect -a Fossil936.Goat.merged.rescaled.bam -b $REGIONS -sorted -g $GENOMEFILE > $OUTDIR/Fossil936.Goat.merged.rescaled.rm.noX.bam
echo "Done subsetting bam file."
date

# Index bam file
echo "Indexing bam file..."
samtools index -@ 2 $OUTDIR/Fossil936.Goat.merged.rescaled.rm.noX.bam
echo "Done indexing bam file."

# Calculate coverage
date
echo "Calculating depth and number of reads at all sites within non-rmMasked regions..."
samtools view -@ 2 -c $OUTDIR/Fossil936.Goat.merged.rescaled.rm.noX.bam | awk '{print "# combined reads " $1}' > $OUTDIR/Fossil936.Goat.merged.rescaled.rm.noX.mapping.results.txt
samtools depth -@ 2 -a -b $REGIONS $OUTDIR/Fossil936.Goat.merged.rescaled.rm.noX.bam | awk '{sum+=$3;cnt++}END{print " Coverage " sum/cnt " Total mapped bp " sum}' >> $OUTDIR/Fossil936.Goat.merged.rescaled.rm.noX.mapping.results.txt
echo "Done calculating depth and number of reads at all sites within non-rmMasked regions."
date

# Picard stats (CollectWgsMetrics)

## Generate Picard intervals list
### Make sequence dictionary for the reference genome
cd /projects/lorenzen/people/KUID/wildebeest/refgenomes
gatk CreateSequenceDictionary -R GCF_001704415.2_ARS1.2_genomic.fasta
### Make intervals from bed file
cd /projects/lorenzen/people/KUID/wildebeest/refgenomes/repeatmasker
gatk BedToIntervalList -I GCF_001704415.2_ARS1.2_genomic_rm.sorted.complement.noX.bed -O GCF_001704415.2_ARS1.2_genomic_rm.sorted.complement.noX.picard.intervals -SD ../GCF_001704415.2_ARS1.2_genomic.dict

## Get coverage stats
INTERVALS="/projects/lorenzen/people/KUID/wildebeest/refgenomes/repeatmasker/GCF_001704415.2_ARS1.2_genomic_rm.sorted.complement.noX.picard.intervals"
REF="/projects/lorenzen/people/KUID/wildebeest/refgenomes/GCF_001704415.2_ARS1.2_genomic.fasta"

echo "Calculating depth within non-rmMasked regions..."
gatk --java-options "-Xmx28G" CollectWgsMetrics -I $OUTDIR/Fossil936.Goat.merged.rescaled.rm.noX.bam \
 -O $OUTDIR/Fossil936.Goat.merged.rescaled.rm.noX.WgsMetrics.txt \
 -R $REF \
 -CAP 1000 \
 --INTERVALS $INTERVALS \
 --READ_LENGTH 41 \
 --COUNT_UNPAIRED \
 --INCLUDE_BQ_HISTOGRAM

date
echo "Done calculating depth within non-rmMasked regions."
