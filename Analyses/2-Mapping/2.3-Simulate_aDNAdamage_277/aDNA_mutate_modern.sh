#!/bin/bash
# Author: Michael V. Westbury
source ~/.bashrc

# Navigate to directory with reads
cd /projects/lorenzen/data/user_owned_folders/Deon_sharing/wildebeest/raw_reads

# Step 1: remove adapters and trim reads to be max 40 and min 30 bp long
~/Software/fastp-0.23.2/bin/fastp --trim_poly_g -i SRR27955277_1.fq.gz --out1 SRR27955277-trimmed.fastq -w 5 --verbose --qualified_quality_phred 25 --length_required 30 --max_len1 40

# Step 2: Add damage patterns
## Uses the tapas sub-programs "filter_fastq" and "multiple_mutate"
conda deactivate
conda activate tapas
~/My_Software/tapas/scripts/filter_fastq --nucleotide @ ~/My_Software/tapas/scripts/multiple_mutate --seed 123 Mutationprofile.tab @ < SRR27955277-trimmed.fastq  > SRR27955277-trimmed_mut.fastq 

# tapas manual: https://mlell.github.io/tapas/index.html
# tapas github: https://github.com/mlell/tapas