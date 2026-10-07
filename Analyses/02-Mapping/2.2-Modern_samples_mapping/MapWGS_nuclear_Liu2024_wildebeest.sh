#!/bin/bash
# Map modern WGS reads

## Modules are loaded in the slurm script used to execute this script
### Modules required: bwa,samtools,gatk4,java-jdk/8.0.112

## Reads should be cleaned before mapping with fastp

## $1 - Threads
## $2 - Clean data path
## $3 - Clean reads parent name (without _1.clean.fq.gz and _2.clean.fq.gz)
## $4 - Code name for sample
## $5 - Results folder
## $6 - Path to reference (including file name)
## $7 - Reference code name

## Sanity check for arguments
echo "Sanity check for provided arguments:"
echo "Threads: $1"
echo "Clean data path: $2"
echo "Clean reads parent name: $3"
echo "Code name for sample: $4"
echo "Results folder: $5"
echo "Reference genome: $6"
echo "Reference genome code name: $7"

cd $5

## Get read group information
## For SRA fastq files, the read group is just going to be the accession ID (SRRXXXXXXX)
### Need: ID,PU,SM,PL,LB (see: https://gatk.broadinstitute.org/hc/en-us/articles/360035890671)
echo "Read Group @RG\tID:$3\tPU:$3\tSM:$4\tPL:ILLUMINA\tLB:$4"


## Index reference genome and map clean reads - produce sorted bam with mate score tags
#bwa index $6
#samtools faidx $6
date "+%F %T"; echo "Mapping reads to reference genome..."
bwa mem -M -a -t $1 \
 -R $(echo "@RG\tID:$3\tPU:$3\tSM:$4\tPL:ILLUMINA\tLB:$4") \
 $6 $2/$4/$3_1.clean.fq.gz $2/$4/$3_2.clean.fq.gz | \
 samtools fixmate -@ $1 -m -u - - | \
 samtools sort -@ $1 -m 7G -O bam -o $5/${4}_${7}.sorted.bam -
date "+%F %T"; echo "Done mapping reads and merging temp bam files."
# Index
date "+%F %T"; echo "Indexing sorted bam file..."
samtools index -@ $1 $5/${4}_${7}.sorted.bam
date "+%F %T"; echo "Done indexing bam file."

## Mark duplicates
date "+%F %T"; echo "Marking duplicate reads..."
#samtools markdup -@ $1 -s -f $5/${4}_${7}.MD.stats.txt -d 2500 $5/${4}_${7}.sorted.bam $5/${4}_${7}.MD.bam
gatk --java-options "-Xmx28G" MarkDuplicates \
 -I $5/${4}_${7}.sorted.bam \
 -O $5/${4}_${7}.MD.bam \
 -M $5/${4}_${7}.MD.metrics.txt \
 --TAGGING_POLICY All \
 --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 # 2500 for NovaSeq reads (uses a patterned flowcell).
MDBAM="$5/${4}_${7}.MD.bam"
date "+%F %T"; echo "Done marking duplicates."

## Filter bam file for high quality mapped reads
date "+%F %T"; echo "Filtering bam file; discarding reads with MQ<30..."
samtools view -@ $1 -q 30 -b --write-index -o $5/${4}_${7}.MD.MQ30.bam $MDBAM
MQ30BAM="$5/${4}_${7}.MD.MQ30.bam"
date "+%F %T"; echo "Done filtering bam file."

## Calculate stats
date "+%F %T"; echo "Calculating mapping statistics with samtools view and depth..."
samtools view -@ $1 -c $MQ30BAM | awk '{print "# combined reads " $1}' > ${MQ30BAM/%.bam/.mapping.results.txt}
samtools depth -a $MQ30BAM | awk '{sum+=$3;cnt++}END{print " Coverage " sum/cnt " Total mapped bp " sum}' >> ${MQ30BAM/%.bam/.mapping.results.txt}
date "+%F %T"; echo "Calculating mapping statistics with gatk CollectWgsMetrics..."
gatk --java-options "-Xmx28G" CollectWgsMetrics \
 -I $MQ30BAM \
 -O ${MQ30BAM/%.bam/.WgsMetrics.txt} \
 -R $6
date "+%F %T"; echo "Done calculating mapping statistics. Remember to remove intermediate bam files."

## Remove intermediate files
rm $5/${4}_${7}.sorted.bam $5/${4}_${7}.sorted.bam.bai
