# Test ANGSD heterozygosity parameters
Compare heterozygosity estimates obtained with different parameters in ANGSD for three versions of the same modern sample (SRR27955277):
 - `277`: High-coverage (~20X) modern; "ground truth"
 - `277_3X`: Low-coverage (~3X) modern; downsampled version of `277`
 - `277_damaged`: Low-coverage (~2.8X) modern sample with ancient DNA damage simulated, based on the Fossil936 damage patterns (see [2.3-Simulate_aDNAdamage_277](https://github.com/DeondeJager/Fincha-Habera-Ancient-Wildebeest/tree/main/Analyses/2-Mapping/2.3-Simulate_aDNAdamage_277))

## Estimate heterozygosity across HiC_scaffold_1 with different parameters
- `angsd_het_SRR27955277_BWD_HiC_scaffold_1_tests.slurm`: Contains the parameters tested for both `277` and `277_3X`
- `angsd_het_SRR27955277_trim_Mut_BWD_HiC_scaffold_1_tests.slurm`: Contains the parameters tested for `277_damaged`
- The output produced from these tests were collated and plotted in R as fig. S6 (see [Figures](https://github.com/DeondeJager/Fincha-Habera-Ancient-Wildebeest/tree/main/Figures))
- This indicated that using `-baq 1` in `ANGSD` was the only parameter that did not produce highly inflated heterozygosity estimates in `277_damaged`.
  - Note that transitions were *always* removed (`-noTrans 1`).

## Statistical comparison of genome-wide heterozygosity
Perform a formal statistical test of heterozygosiy across the entire genome for the three versions of the sample, using the `-baq 1` parameter and with transitions removed (as always) `noTrans 1`.
- `angsd_het_SRR27955277_BWD_HiC.slurm` = `277`
- `angsd_het_SRR27955277_3X_BWD_HiC.slurm` = `277_3X`
- `angsd_het_SRR27955277_trim_Mut_BWD_HiC.slurm` = `277_damaged`

The output produced from these test were collated and plotted in R as fig. S7, which is where the Wilcoxon signed-rank test for a difference in means was also performed (see Figures](https://github.com/DeondeJager/Fincha-Habera-Ancient-Wildebeest/tree/main/Figures))
