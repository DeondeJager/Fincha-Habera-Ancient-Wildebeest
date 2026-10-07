#!/bin/bash

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules and set up enviroment variables
module purge
module load getorganelle/1.7.7.0
READS=/projects/lorenzen/people/KUID/wildebeest/clean_data/Liu_2024
SEED=/projects/lorenzen/people/KUID/wildebeest/refgenomes
OUTDIR=/projects/lorenzen/people/KUID/wildebeest/GetOrganelle

# Navigate to working directory
cd $OUTDIR

## Without seed
echo "Assembling SRR27955239 mitogenome WITHOUT seed..."
get_organelle_from_reads.py -1 $READS/SRR27955239/SRR27955239_1.clean.fq.gz -2 $READS/SRR27955239/SRR27955239_2.clean.fq.gz -F animal_mt -o $OUTDIR/SRR27955239 -R 10 -t 4
echo "Done assembling SRR27955239 mitogenome WITHOUT seed."
date


