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
OUTDIR=/projects/lorenzen/people/KUID/wildebeest/admixture_graph/BWD_HiC

# Generate a pseudo-haploid file
## n81_modern is without the fossil, and only includes the modern samples including SRR27955277_trim_Mut but excluding SRR27955277, also has the outgroups (hartebeest)
## Note: Transitions removed
cd $OUTDIR
angsd -b $OUTDIR/bamlist_n81_modern_BWD_HiC.txt \
 -dohaplocall 1 -minMinor 4 -maxMis 9 \
 -doCounts 1 -minMapQ 30 -minQ 20 -baq 1 -noTrans 1 \
 -r HiC_scaffold_24: \
 -remove_bads 1 -uniqueOnly 1 \
 -setMinDepth 461 -setMaxDepth 2854 \
 -nThreads 20 \
 -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta \
 -out $OUTDIR/n81_modern_HiC_scaffold_24

