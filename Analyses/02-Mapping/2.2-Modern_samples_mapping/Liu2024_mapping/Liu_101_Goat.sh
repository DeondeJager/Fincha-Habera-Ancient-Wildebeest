#!/bin/bash

# Map wildebeest SRR27955318 reads

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules (conda environments are loaded and unloaded as need the *.sh script)
module load bwa/0.7.17
module load samtools/1.15
module load gatk4/4.3.0.0

# Navigate to working directory
cd /projects/lorenzen/people/KUID/wildebeest/clean_data/Liu_2024/

# Order of input arguments
## $1 - Threads
## $2 - Clean data path
## $3 - Clean reads parent name (without _1.clean.fq.gz and _2.clean.fq.gz)
## $4 - Code name for sample
## $5 - Results folder
## $6 - Path to reference (including file name)
## $7 - Reference code name

# Call script
## Map with bwa
bash /home/KUID/scripts/mapping/MapWGS_nuclear_Liu2024_wildebeest.sh \
 10 \
 /projects/lorenzen/people/KUID/wildebeest/clean_data/Liu_2024 \
 SRR27955318 \
 SRR27955318 \
 /projects/lorenzen/people/KUID/wildebeest/bams \
 /projects/lorenzen/people/KUID/wildebeest/refgenomes/GCF_001704415.2_ARS1.2_genomic.fasta \
 Goat

