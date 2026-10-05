#!/bin/bash

# Subset bam file of sample SRR27955240 to only goat autosomes (no X) and with repeat regions masked

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load module
module purge
module load bedtools/2.31.0
module load samtools/1.15
module load java-jdk/8.0.112 # For gatk
module load gatk4/4.3.0.0 # For CollectWgsMetrics

#--- Goat reference genome ---#
REGIONS="/projects/lorenzen/people/KUID/wildebeest/refgenomes/repeatmasker/GCF_001704415.2_ARS1.2_genomic_rm.sorted.complement.noX.bed"
GENOMEFILE="/projects/lorenzen/people/KUID/wildebeest/refgenomes/GCF_001704415.2_ARS1.2_genomic.chr.sizes.txt"
OUTDIR=/projects/lorenzen/people/KUID/wildebeest/bams

# Navigate to working directory
cd $OUTDIR

# Subset bam file using regions file
date
echo "Subsetting bam file..."
bedtools intersect -a SRR27955240_Goat.MD.MQ30.bam -b $REGIONS -sorted -g $GENOMEFILE > $OUTDIR/SRR27955240_Goat.MD.MQ30.rm.noX.bam
echo "Done subsetting bam file."
date

# Index bam file
cd $OUTDIR
echo "Indexing bam file..."
samtools index -@ 2 SRR27955240_Goat.MD.MQ30.rm.noX.bam
echo "Done indexing bam file."

# Calculate coverage
date
echo "Calculating depth and number of reads at all sites within non-rmMasked regions..."
samtools view -@ 2 -c SRR27955240_Goat.MD.MQ30.rm.noX.bam | awk '{print "# combined reads " $1}' > SRR27955240_Goat.MD.MQ30.rm.noX.mapping.results.txt
samtools depth -@ 2 -a -b $REGIONS SRR27955240_Goat.MD.MQ30.rm.noX.bam | awk '{sum+=$3;cnt++}END{print " Coverage " sum/cnt " Total mapped bp " sum}' >> SRR27955240_Goat.MD.MQ30.rm.noX.mapping.results.txt
echo "Done calculating depth and number of reads at all sites within non-rmMasked regions."
date

# Picard stats (CollectWgsMetrics)
INTERVALS="/projects/lorenzen/people/KUID/wildebeest/refgenomes/repeatmasker/GCF_001704415.2_ARS1.2_genomic_rm.sorted.complement.noX.picard.intervals"
REF="/projects/lorenzen/people/KUID/wildebeest/refgenomes/GCF_001704415.2_ARS1.2_genomic.fasta"

echo "Calculating depth within non-rmMasked regions..."
gatk --java-options "-Xmx28G" CollectWgsMetrics -I SRR27955240_Goat.MD.MQ30.rm.noX.bam \
 -O SRR27955240_Goat.MD.MQ30.rm.noX.WgsMetrics.txt \
 -R $REF \
 -CAP 1000 \
 --INTERVALS $INTERVALS \
 --INCLUDE_BQ_HISTOGRAM

date
echo "Done calculating depth within non-rmMasked regions."
