#!/bin/bash

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used 
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules
module load angsd/0.921

# Navigate to working directory and set up variables
OUT=/projects/lorenzen/people/KUID/wildebeest/het/BWD_HiC
BAMS=/projects/lorenzen/people/KUID/wildebeest/bams
REF=/projects/lorenzen/people/KUID/wildebeest/refgenomes
ANC=/projects/lorenzen/people/KUID/wildebeest/fasta

# Insert commands here:
cd $OUT
## Create saf.idx file
### Max depth set to: Median_cov + 2*SD_cov (as calculated by Picard's CollectWgsMetrics) - different for each sample
### Note: Using filtered bam files with mapping quality of >=30, but mapping qual is recalibrates, so -minMapQ30 is added
## No transitions, unfolded sfs, as ancenstral state is given by hartebeest fasta
angsd -i $BAMS/SRR27955277_BWD_HiC-Human.MD.MQ30.rm.noX.bam -rf $REF/BWD_HiC.chr.fasta.autosomes.ANGSD.txt -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 1 -doCounts 1 -setMinDepth 7 -setMaxDepth 40 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -nThreads 10 -out $OUT/SRR27955277_BWD_HiC_noTrans
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_BWD_HiC_noTrans.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_BWD_HiC_noTrans.10Mb.ml
