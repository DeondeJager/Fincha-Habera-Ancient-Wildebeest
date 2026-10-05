#!/bin/bash
# Make BED regions file of repeat (low complexity) regions to mask in bam files
## Author: Deon de Jager. 20 March 2024.


## Steps overview:
### 1. Run repeatmasker on reference genome
### 2. Convert repeatmasker output to BED file
### 3. Complement this BED file to get the parts of the genome we want

## Steps
### 1. See repeatmasker script (completed before).

### 2. Convert repeatmasker output to BED file
#### Using the rmsk2bed function from bedops and keep only the first three columns (cut)
#### https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/rmsk2bed.html
#### Start an interactive session
cd /projects/lorenzen/people/KUID/wildebeest/refgenomes/repeatmasker
module load bedops/2.4.41
rmsk2bed < GCF_001704415.2_ARS1.2_genomic.fasta.out | cut -f1,2,3 > GCF_001704415.2_ARS1.2_genomic_rm.bed

### 3. Sort and complement this BED file to get the parts of the genome we want
#### Here, we use bedtools, since its complement function is slightly easier to use than bedops'
#### https://bedtools.readthedocs.io/en/latest/content/tools/complement.html
#### First, we need a file that lists the chromosome names and their lengths,
#### which we can get by indexing the ref genome using samtools faidx:
cd ../
module load samtools/1.15
samtools faidx GCF_001704415.2_ARS1.2_genomic.fasta
#### And then extracting the first two columns:
cut -f1,2 GCF_001704415.2_ARS1.2_genomic.fasta.fai > GCF_001704415.2_ARS1.2_genomic.chr.sizes.txt
#### Then sort the bed file to be in the same order as the chromosomes in the genome (required later by angsd)
module load bedtools/2.31.0
cd repeatmasker/
bedtools sort -faidx ../GCF_001704415.2_ARS1.2_genomic.fasta.fai -i GCF_001704415.2_ARS1.2_genomic_rm.bed > GCF_001704415.2_ARS1.2_genomic_rm.sorted.bed
#### Then complement the repeatmasker BED file
bedtools complement -i GCF_001704415.2_ARS1.2_genomic_rm.sorted.bed -g ../GCF_001704415.2_ARS1.2_genomic.chr.sizes.txt > GCF_001704415.2_ARS1.2_genomic_rm.sorted.complement.bed

### 4. Remove X chromosome scaffolds from the bed file
#### The X chromosome is represented by two large scaffods: NW_017189516.1 and NW_017189517.1
grep -v -F 'NW_017189516.1' GCF_001704415.2_ARS1.2_genomic_rm.sorted.complement.bed | grep -v -F 'NW_017189517.1' > GCF_001704415.2_ARS1.2_genomic_rm.sorted.complement.noX.bed
# This file is then used in the subset_BAM_byRegions_Goat.slurm script
