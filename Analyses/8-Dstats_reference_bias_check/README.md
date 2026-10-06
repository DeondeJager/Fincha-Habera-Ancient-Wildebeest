# D-statistics (ABBA-BABA test) for reference bias
To investigate potential reference genome bias, we used the simulated sample, `277_damaged`, and conducted an ABBA-BABA test in ANGSD with this sample and its undamaged, high-coverage counterpart, `277`.

This was done for the BWD_HiC and goat reference genomes.

## BWD_HiC reference
- `angsd_Abbababba_BWD_HiC_refBiasCheck.slurm`
  - Uses both the `bamlist_n10_BWD_HiC_refBiasCheck.txt` and `indNames_noOutgroup_refBiasCheck.txt` files.
  
## Goat reference
- `angsd_Abbababba_Goat_refBiasCheck.slurm`
  - Uses both the `bamlist_n10_Goat_refBiasCheck.txt` and `indNames_noOutgroup_refBiasCheck.txt` files.

## Plotting
The `*.jackknife.txt` output files from the above analyses were then used to generate the D-statistics plots in fig. S9 (see [Figures](https://github.com/DeondeJager/Fincha-Habera-Ancient-Wildebeest/tree/main/Figures)), which indicated the presence of reference bias when the goat genome was used as the reference, but not when the wildebeest BWD_HiC genome was used.
