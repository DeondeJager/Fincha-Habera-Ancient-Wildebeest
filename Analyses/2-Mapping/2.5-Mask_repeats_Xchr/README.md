# Mask repeat regions and remove chromosome X of mapDamage-rescaled bam files

## Wildebeest BWD_HiC bams
1. Make a `BED` file with regions to keep (the rest will be masked): See `prep_RepeatRegions_BWD_HiC.sh`.
2. Retain only regions _not_ in the `BWD_HiC.repeatmasker.trf.windowmasker.chr.sort.complement.noX.bed` file: See `subset_BAM_byRegions_FossilWildebeest.slurm`

## Goat bams
1. Run `repeatmasker` on goat reference genome to identify repeat regions: See `repeatmasker_goat.slurm`
2. Make `BED` file with regions to keep (the rest will be masked): See `prep_RepeatRegions_Goat.sh`
3. Retain only regions _not_ in the `GCF_001704415.2_ARS1.2_genomic_rm.sorted.complement.noX.bed` file: See `subset_BAM_byRegions_Goat.slurm`

## Modern wildebeest samples
The bam files of the modern wildebeest samples from Liu et al. (2024) were subset in the same way.
 - See the subfolder `Liu2024_mask` for the scripts with the actual commands.
 - The `Liu2024_wildebeest.batch.slurm` and `Liu2024_wildebeest.Goat.batch.slurm` slurm scripts submits these as array jobs on the cluster. 
