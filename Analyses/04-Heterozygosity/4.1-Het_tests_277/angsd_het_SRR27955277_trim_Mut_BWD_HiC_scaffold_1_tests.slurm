#!/bin/bash
#SBATCH --job-name=het_SRR27955277_trim_Mut_BWD_HiC_tests     # Job name
#SBATCH --ntasks=1                     # Number of tasks per node
#SBATCH --cpus-per-task=10              # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G              # Memory required per CPU
#SBATCH --time=04:00:00                # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/het/BWD_HiC/tests/angsd_het_SRR27955277_trim_Mut_BWD_HiC_tests_%j.out    # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/het/BWD_HiC/tests/angsd_het_SRR27955277_trim_Mut_BWD_HiC_tests_%j.err     # Standard error log

# The %j in the --output & --error lines tells SLURM to substitute the job ID in the name of the output file
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used 
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules
module load angsd/0.921

# Navigate to working directory and set up variables
OUT=/projects/lorenzen/people/KUID/wildebeest/het/BWD_HiC/tests
BAMS=/projects/lorenzen/people/KUID/wildebeest/paleomix/competitive_mapping/SRR27955277/BWD_HiC/mapDamage
REF=/projects/lorenzen/people/KUID/wildebeest/refgenomes
ANC=/projects/lorenzen/people/KUID/wildebeest/fasta

# Insert commands here:
cd $OUT
## Create saf.idx file
### Max depth set to: Median_cov + 2*SD_cov (as calculated by Picard's CollectWgsMetrics) - different for each sample
### Note: Using filtered bam files with mapping quality of >=30, but mapping qual is recalibrated, so -minMapQ30 is added

#--- 1. Initial estimation - same params as modern sample (baq2) ---#
# Result was much higher het than non-damaged version of this sample, so need to adjust parameters until I get a similar estimate, and then use those parameters for the fossil.
## No transitions, unfolded sfs, as ancenstral state is given by hartebeest fasta
angsd -i $BAMS/SRR27955277_trim_Mut.BWD_HiC-Human.rescaled.rm.noX.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 2 -doCounts 1 -setMinDepth 1 -setMaxDepth 7 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -nThreads 10 -out $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq2
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq2.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq2.10Mb.ml
## Het is much higher than the undamaged modern sample with baq2. It appears there is something not quite right with baq2 and ancient samples.

#--- 2. Discard (trim) 10 bases at the ends of reads ---#
# Do super conservative trim of 10 bases initially (based on mapDamage plots) to see if it has an effect. 
# If it has an effect, then do trims from 1-9 to find a balance between reducing error and retaining information.
## Restrict testing to longest scaffold (HiC_scaffold_1) for faster results.
angsd -i $BAMS/SRR27955277_trim_Mut.BWD_HiC-Human.rescaled.rm.noX.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -trim 10 -doSaf 1 -noTrans 1 -GL 2 -baq 2 -doCounts 1 -setMinDepth 1 -setMaxDepth 7 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -nThreads 10 -out $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_trim10
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_trim10.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_trim10.10Mb.ml
## Het is slightly higher than in 1.

#--- 3. Adjust mapQ for excessive mismatches (-C 50) ---#
angsd -i $BAMS/SRR27955277_trim_Mut.BWD_HiC-Human.rescaled.rm.noX.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 2 -C 50 -doCounts 1 -setMinDepth 1 -setMaxDepth 7 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -nThreads 10 -out $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_C50
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_C50.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_C50.10Mb.ml
## Het is about the same as in 1.

#--- 4. Discard (trim) 4 bases at the ends of reads ---#
# trim10 increased het. So let's try a less conservative trim and see what happens.
angsd -i $BAMS/SRR27955277_trim_Mut.BWD_HiC-Human.rescaled.rm.noX.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -trim 4 -doSaf 1 -noTrans 1 -GL 2 -baq 2 -doCounts 1 -setMinDepth 1 -setMaxDepth 7 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -nThreads 10 -out $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_trim4
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_trim4.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_trim4.10Mb.ml
## Het is slighlty lower than trim10, but still higher than in 1.

#--- 5. Use baq 1 ---#
angsd -i $BAMS/SRR27955277_trim_Mut.BWD_HiC-Human.rescaled.rm.noX.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 1 -doCounts 1 -setMinDepth 1 -setMaxDepth 7 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -nThreads 10 -out $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq1
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq1.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq1.10Mb.ml
## Het is much lower than in 1 and equivalent to the undamaged modern sample with baq1. The 3X modern sample with baq1 is slightly lower than both of the aformentioned. 

#--- 6. Use folded spectrum ---#
# The sfs is folded by realSFS.
angsd -i $BAMS/SRR27955277_trim_Mut.BWD_HiC-Human.rescaled.rm.noX.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $REF/BWD_HiC.chr_Human_fossilMito.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 2 -doCounts 1 -setMinDepth 1 -setMaxDepth 7 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -nThreads 10 -out $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_fold
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_fold.saf.idx -fold 1 -nSites 10000000 -P 10 > $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_fold.10Mb.ml
## Het is about the same as in 1.

#--- 7. No baq ---#
angsd -i $BAMS/SRR27955277_trim_Mut.BWD_HiC-Human.rescaled.rm.noX.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 0 -doCounts 1 -setMinDepth 1 -setMaxDepth 7 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -nThreads 10 -out $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq0
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq0.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq0.10Mb.ml

#--- 8. baq 2 & minDepth 3 ---#
angsd -i $BAMS/SRR27955277_trim_Mut.BWD_HiC-Human.rescaled.rm.noX.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 2 -doCounts 1 -setMinDepth 3 -setMaxDepth 7 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -nThreads 10 -out $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq2_minD3
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq2_minD3.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq2_minD3.10Mb.ml

#--- 9. baq 1 & minDepth 3 ---#
angsd -i $BAMS/SRR27955277_trim_Mut.BWD_HiC-Human.rescaled.rm.noX.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 1 -doCounts 1 -setMinDepth 3 -setMaxDepth 7 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -nThreads 10 -out $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq1_minD3
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq1_minD3.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq1_minD3.10Mb.ml

#--- 10. baq 0 & minDepth 3 ---#
angsd -i $BAMS/SRR27955277_trim_Mut.BWD_HiC-Human.rescaled.rm.noX.bam -r HiC_scaffold_1: -ref $REF/BWD_HiC.chr_Human_fossilMito.fasta -anc $ANC/SRR27955304_BWD_HiC-Human.MD.MQ30.rm.noX.fasta -doSaf 1 -noTrans 1 -GL 2 -baq 0 -doCounts 1 -setMinDepth 3 -setMaxDepth 7 -minQ 20 -minMapQ 30 -uniqueOnly 1 -remove_bads 1 -nThreads 10 -out $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq0_minD3
## This is followed by the actual estimation using the realSFS subprogram of ANGSD in 10 Mb windows
realSFS $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq0_minD3.saf.idx -nSites 10000000 -P 10 > $OUT/SRR27955277_trim_Mut_BWD_HiC_noTrans_baq0_minD3.10Mb.ml

