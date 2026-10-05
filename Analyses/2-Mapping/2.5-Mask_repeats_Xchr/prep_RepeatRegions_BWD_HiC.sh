#!/bin/bash
# Make BED regions file of repeat (low complexity) regions to mask in bam files
## Author: Deon de Jager. 12 September 2024.

## Steps overview:
### 1. Obtain repeat file for BWD_HiC assembly
### 2. Convert repeat file (gff format) to BED file
### 3. Complement this BED file to get the parts of the genome we want
### 4. Remove the HiC scaffold that corresponds to the X chromosome

## Steps
### 1. Download the repeatmasker & windowmasker file from here: https://dnazoo.s3.wasabisys.com/index.html?prefix=Connochaetes_taurinus/
####	Filename: BWD_HiC.repeatmasker.trf.windowmasker.gff.gz. Put it in the repeatmasker directory
cd /projects/lorenzen/people/KUID/wildebeest/refgenomes/repeatmasker
gzip -d BWD_HiC.repeatmasker.trf.windowmasker.gff.gz # Unzip

### 2. Convert repeat gff file to BED file
####	2.1. Convert gff to bed using bedops
module load bedops/2.4.41
gff2bed < BWD_HiC.repeatmasker.trf.windowmasker.gff > BWD_HiC.repeatmasker.trf.windowmasker.bed
####	2.2. Get a list of the BWD_HiC scaffolds from the .fai index of the reference genome and name it BWD_HiC.chr.txt
module load samtools/1.15
samtools faidx ../BWD_HiC.chr.fasta
cut -f1 ../BWD_HiC.chr.fasta.fai > BWD_HiC.chr.txt
####	2.3. Filter repeat bed file to only contain lines corresponding to the chromosome scaffolds in BWD_HiC.chr.txt
grep -f BWD_HiC.chr.txt -w BWD_HiC.repeatmasker.trf.windowmasker.bed | awk '{print $1"\t"$2"\t"$3}' > BWD_HiC.repeatmasker.trf.windowmasker.chr.bed
			# -f: get pattern to grep from file
			# -w: match whole words only (otherwise other scaffolds like "HiC_scaffold_111" is also included for the pattern "HiC_scaffold_1"
			# awk: This part prints only columns 1, 2, and 3 from the bed file, which is the scaffold name, start, and end positions.
####	2.4. Sort the bed file to match the order of scaffolds in the ref genome, as the current order is HiC_scaffold_1, HiC_scaffold_10, etc. instead of HiC_scaffold_1, HiC_scaffold_2, etc. You can check the order as follows:
cut -f1 BWD_HiC.repeatmasker.trf.windowmasker.chr.bed | uniq
sort -V BWD_HiC.repeatmasker.trf.windowmasker.chr.bed > BWD_HiC.repeatmasker.trf.windowmasker.chr.sort.bed
			# -V: natural sort of (version) numbers within text.
cut -f1 BWD_HiC.repeatmasker.trf.windowmasker.chr.sort.bed | uniq # Check the sorting worked as expected

### 3. Complement this BED file to get the parts of the genome we want
#### Here, we use bedtools, since its complement function is slightly easier to use than bedops'
#### https://bedtools.readthedocs.io/en/latest/content/tools/complement.html
#### 	3.1. First, we need a "genome" file that lists the chromosome names and their lengths, which we can get from the .fai file generated with samtools in step 2.2 by extracting the first two columns:
cut -f1,2 ../BWD_HiC.chr.fasta.fai > ../BWD_HiC.chr.fasta.chr.sizes.txt
####	3.2. Then complement the bed file using bedtools and the "genome" file generated in step 3.1
module load bedtools/2.31.0
bedtools complement -i BWD_HiC.repeatmasker.trf.windowmasker.chr.sort.bed -g ../BWD_HiC.chr.fasta.chr.sizes.txt > BWD_HiC.repeatmasker.trf.windowmasker.chr.sort.complement.bed
####	3.3. You can check the lengths of each bed file (the one indicating the repeat regions, and the complement indicating non-repeat regions)
cat BWD_HiC.repeatmasker.trf.windowmasker.chr.sort.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' # Output: 2,079,380,334 bp in repeats
cat BWD_HiC.repeatmasker.trf.windowmasker.chr.sort.complement.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' # 1,189,713,602 bp not in repeats
# Thus, a larger part of the genome is repeats than is not repeats. I think this is expected/relatively normal.

####	3.4. Remove intermediate bed files
rm BWD_HiC.repeatmasker.trf.windowmasker.bed BWD_HiC.repeatmasker.trf.windowmasker.chr.bed

### 4. Remove the HiC scaffold that corresponds to the X chromosome
#### This was done by a combination of three approaches:
#### 4.1. Searching the gff file from DNA zoo for AMEL (for amelogenin, an X- and Y-chr gene):
####		- This gave the result of matches for AMELX and AMELY (and others) on HiC_scaffold_3 and HiC_scaffold_22
#### 4.2. Running SatsumaSynteny using the bovine X and Y chromosomes (see the directory ~/data/wildebeest/SeXY/BWD_HiC)
####		- This gave the result of mainly HiC_scaffold_3 for both the X and Y chr, but also some other regions scattered across other HiC_scaffolds. See files X.bed and Y.bed in the SeXY/BWD_HiC directory.
####		- Notably, HiC_scaffold_22 was absent from both the X and Y SatsumaSynteny results. See files X.bed.scaffolds.counts and Y.bed.scaffolds.counts
#### 4.3. Coverage was calculated per scaffold. This was done per-site and then the mean was calculated for 10M random sites 10 times for each scaffold (based on the CoverageCalculation.sh script from SeXY). See script ~/scripts/samtools/fossilWildebeest/CoverageCalculation_BWD-HiC.slurm.
####		- This was done for a male individual from Liu et al. (2024) (SRR27955227) with the idea that the X-chr would have half the coverage of the rest of the chromosomes (and the Y chr, if it is present - the sex of the reference genome individual is unknown at this point). A female individual was used as a control (SRR27955224).
####		- The results were plotted in R and showed that HiC_scaffold_3 had half the coverage of the rest of the scaffolds. See plot male_female_scaffold_coverage.pdf.
#### As a result, it was decided that HiC_scaffold_3 was the (only) one that corresponds to a sex chromosome (the X chr).

### 5. Remove the X-chr HiC_scaffold_3 from the repeatmasker bed file, so that it is excluded when subsetting the bams.
grep -v -F 'HiC_scaffold_3' BWD_HiC.repeatmasker.trf.windowmasker.chr.sort.complement.bed > BWD_HiC.repeatmasker.trf.windowmasker.chr.sort.complement.noX.bed
# This file is then used in the subset_BAM_byRegions_FossilWildebeest.slurm script