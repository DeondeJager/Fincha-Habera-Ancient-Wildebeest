# Generate pseudohaploid datasets
- The scripts in the folders `BWD_HiC` and `Goat` generate the pseudohaploid dataset using `ANGSD` for each chromosome separately.
  - These scripts require the `bamlist_*.txt` files, which gives the path to all the bam files to be included in the pseudohaploid dataset.
- The `.slurm` scripts submit these as slurm arrays and also contain commands to join the per-chromosome `.haplo` files into one.
- There are two datasets as indicated by the filenames:
  - `*n82_fossil*`: Contains Fossil936 and the modern samples (excludes 277_damaged, but includes its undamaged bam, 277) + hartebeest outgroups.
  - `*n81_modern*`: Excludes Fossil936 and includes only modern samples, but with 277_damaged (aka SRR27955277_trim_Mut) and excluding its undamaged counterpart SRR27955277; also has the hartebeest outgroups.
