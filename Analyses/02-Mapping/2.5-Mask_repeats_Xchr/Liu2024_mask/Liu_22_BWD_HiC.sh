#!/bin/bash

# Subset bam file of sample SRR27955233 to only BWD_HiC autosomes (no X/HiC_scaffold_3) and with repeat regions masked

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

#--- BWD_HiC-Human reference genome ---#
REGIONS="/projects/lorenzen/people/KUID/wildebeest/refgenomes/repeatmasker/BWD_HiC.repeatmasker.trf.windowmasker.chr.sort.complement.noX.bed"
GENOMEFILE="/projects/lorenzen/people/KUID/wildebeest/refgenomes/BWD_HiC.chr_Human_fossilMito.fasta.chr.sizes.txt"
OUTDIR=/projects/lorenzen/people/KUID/wildebeest/bams

# Navigate to working directory
cd $OUTDIR

# Subset bam file using regions file
date
echo "Subsetting bam file..."
bedtools intersect -a SRR27955233_BWD_HiC-Human.MD.MQ30.bam -b $REGIONS -sorted -g $GENOMEFILE > $OUTDIR/SRR27955233_BWD_HiC-Human.MD.MQ30.rm.noX.bam
echo "Done subsetting bam file."
date

# Index bam file
cd $OUTDIR
echo "Indexing bam file..."
samtools index -@ 2 SRR27955233_BWD_HiC-Human.MD.MQ30.rm.noX.bam
echo "Done indexing bam file."

# Calculate coverage
date
echo "Calculating depth and number of reads at all sites within non-rmMasked regions..."
samtools view -@ 2 -c SRR27955233_BWD_HiC-Human.MD.MQ30.rm.noX.bam | awk '{print "# combined reads " $1}' > SRR27955233_BWD_HiC-Human.MD.MQ30.rm.noX.mapping.results.txt
samtools depth -@ 2 -a -b $REGIONS SRR27955233_BWD_HiC-Human.MD.MQ30.rm.noX.bam | awk '{sum+=$3;cnt++}END{print " Coverage " sum/cnt " Total mapped bp " sum}' >> SRR27955233_BWD_HiC-Human.MD.MQ30.rm.noX.mapping.results.txt
echo "Done calculating depth and number of reads at all sites within non-rmMasked regions."
date

# Picard stats (CollectWgsMetrics)
## Generate Picard intervals list - see picard_intervals.sh script
INTERVALS="/projects/lorenzen/people/KUID/wildebeest/refgenomes/repeatmasker/BWD_HiC.repeatmasker.trf.windowmasker.chr.sort.complement.noX.intervals.picard"
REF="/projects/lorenzen/people/KUID/wildebeest/refgenomes/BWD_HiC.chr_Human_fossilMito.fasta"

echo "Calculating depth within non-rmMasked regions..."
gatk --java-options "-Xmx28G" CollectWgsMetrics -I SRR27955233_BWD_HiC-Human.MD.MQ30.rm.noX.bam \
 -O SRR27955233_BWD_HiC-Human.MD.MQ30.rm.noX.WgsMetrics.txt \
 -R $REF \
 -CAP 1000 \
 --INTERVALS $INTERVALS \
 --INCLUDE_BQ_HISTOGRAM

date
echo "Done calculating depth within non-rmMasked regions."
