#!/bin/bash

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules
module purge
module load angsd/0.921
REF=/projects/lorenzen/people/KUID/wildebeest/refgenomes
OUTDIR=/projects/lorenzen/people/KUID/wildebeest/admixture_graph/Goat

# Generate a pseudo-haploid file for n81 dataset in angsd (modern homogenous pops+Fossil936)
## Note: Transitions removed
cd $OUTDIR
angsd -b $OUTDIR/bamlist_n81_modern_Goat.txt -dohaplocall 1 \
 -doCounts 1 -minMapQ 30 -minQ 20 -baq 1 \
 -noTrans 1 \
 -r NC_030832.1: \
 -minMinor 4 -maxMis 9 \
 -remove_bads 1 -uniqueOnly 1 \
 -setMinDepth 431 -setMaxDepth 3323 \
 -nThreads 20 \
 -ref $REF/GCF_001704415.2_ARS1.2_genomic.fasta \
 -out $OUTDIR/n81_modern_NC_030832.1

