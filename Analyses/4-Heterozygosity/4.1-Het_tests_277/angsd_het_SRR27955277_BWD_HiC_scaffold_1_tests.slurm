#!/bin/bash
#SBATCH --job-name=het_SRR27955277_BWD_HiC_tests     # Job name
#SBATCH --ntasks=1                     # Number of tasks per node
#SBATCH --cpus-per-task=10              # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G              # Memory required per CPU
#SBATCH --time=10:00:00                # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/het/BWD_HiC/tests/angsd_het_SRR27955277_BWD_HiC_tests_%j.out    # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/het/BWD_HiC/tests/angsd_het_SRR27955277_BWD_HiC_tests_%j.err     # Standard error log

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used 
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules
module load angsd/0.921

# Navigate to working directory and set up variables
OUT=/projects/lorenzen/people/KUID/wildebeest/het/BWD_HiC/tests
BAMS=/projects/lorenzen/people/KUID/wildebeest/bams
REF=/projects/lorenzen/people/KUID/wildebeest/refgenomes
ANC=/projects/lorenzen/people/KUID/wildebeest/fasta

# Insert commands here:
cd $OUT
## Create saf.idx file
### Max depth set to: Median_cov + 2*SD_cov (as calculated by Picard's CollectWgsMetrics) - different for each sample
### Note: Using filtered bam files with mapping quality of >=30, but mapping qual is recalibrates, so -minMapQ30 is added

# Restrict to longest scaffold for time reasons (HiC_scaffold_1)

#--- 1. Original het estimate (baq2) ---#
## No transitions, unfolded sfs, as ancenstral state is given by hartebeest fasta
angsd -i $BAMS/SRR27955277_BWD_HiC-Human.MD.MQ30.rm.noX.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 2 -doCounts 1 -setMinDepth 7 -setMaxDepth 40 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -nThreads 10 -out $OUT/SRR27955277_BWD_HiC_noTrans_baq2
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_BWD_HiC_noTrans_baq2.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_BWD_HiC_noTrans_baq2.10Mb.ml

#--- 2. Original het estimate (baq2) on 3X bam ---#
angsd -i $BAMS/BWD_HiC_3X/SRR27955277_BWD_HiC-Human.MD.MQ30.rm.noX.3X.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 2 -doCounts 1 -setMinDepth 1 -setMaxDepth 8 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -nThreads 10 -out $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq2
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq2.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq2.10Mb.ml

#--- 3. Using baq1 ---#
angsd -i $BAMS/SRR27955277_BWD_HiC-Human.MD.MQ30.rm.noX.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 1 -doCounts 1 -setMinDepth 7 -setMaxDepth 40 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -nThreads 10 -out $OUT/SRR27955277_BWD_HiC_noTrans_baq1
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_BWD_HiC_noTrans_baq1.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_BWD_HiC_noTrans_baq1.10Mb.ml

#--- 4. Using baq1 on 3X bam ---#
angsd -i $BAMS/BWD_HiC_3X/SRR27955277_BWD_HiC-Human.MD.MQ30.rm.noX.3X.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 1 -doCounts 1 -setMinDepth 1 -setMaxDepth 8 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -nThreads 10 -out $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq1
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq1.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq1.10Mb.ml

#--- 5. No baq high depth bam ---#
angsd -i $BAMS/SRR27955277_BWD_HiC-Human.MD.MQ30.rm.noX.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 0 -doCounts 1 -setMinDepth 7 -setMaxDepth 40 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -nThreads 10 -out $OUT/SRR27955277_BWD_HiC_noTrans_baq0
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_BWD_HiC_noTrans_baq0.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_BWD_HiC_noTrans_baq0.10Mb.ml

#--- 6. No baq 3X bam ---#
angsd -i $BAMS/BWD_HiC_3X/SRR27955277_BWD_HiC-Human.MD.MQ30.rm.noX.3X.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 0 -doCounts 1 -setMinDepth 1 -setMaxDepth 8 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -nThreads 10 -out $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq0
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq0.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq0.10Mb.ml

#--- 7. baq 2 with minDepth 3 on 3X bam ---#
angsd -i $BAMS/BWD_HiC_3X/SRR27955277_BWD_HiC-Human.MD.MQ30.rm.noX.3X.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 2 -doCounts 1 -setMinDepth 3 -setMaxDepth 8 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -nThreads 10 -out $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq2_minD3
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq2_minD3.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq2_minD3.10Mb.ml

#--- 8. baq 1 with minDepth 3 on 3X bam ---#
angsd -i $BAMS/BWD_HiC_3X/SRR27955277_BWD_HiC-Human.MD.MQ30.rm.noX.3X.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 1 -doCounts 1 -setMinDepth 3 -setMaxDepth 8 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -nThreads 10 -out $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq1_minD3
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq1_minD3.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq1_minD3.10Mb.ml

#--- 9. baq 0 with minDepth 3 on 3X bam ---#
angsd -i $BAMS/BWD_HiC_3X/SRR27955277_BWD_HiC-Human.MD.MQ30.rm.noX.3X.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 0 -doCounts 1 -setMinDepth 3 -setMaxDepth 8 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -nThreads 10 -out $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq0_minD3
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq0_minD3.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_BWD_HiC_noTrans_3X_baq0_minD3.10Mb.ml

