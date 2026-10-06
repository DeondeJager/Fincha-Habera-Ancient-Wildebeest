# Estimate heterozygosity for Fossil936 and contemporary samples
Based on the results of the heterozygosity testing done in `4.1-Het_tests_277`, we estimated genome-wide heterozygosity for all samples using `ANGSD` with the `baq 1` option, as well as no transitions `-noTrans 1`, amongst other filters (see scripts).

## Fossil936
- `angsd_het_Fossil936_BWD_HiC.slurm`

## Contemporary samples
See folder `Liu2024_het` for per-sample commands.
The scripts are submitted as slurm array jobs by `angsd_het_Liu2024_BWD_HiC.batch.slurm`.

Of note `-setMinDepth` and `-setMaxDepth` were individually set for each sample based on the median and standard deviation of its coverage, as calculated by Picard’s CollectWgsMetrics in the mapping script ([2.2-Modern_samples_mapping/MapWGS_nuclear_Liu2024_wildebeest.sh](https://github.com/DeondeJager/Fincha-Habera-Ancient-Wildebeest/tree/main/Analyses/2-Mapping/2.2-Modern_samples_mapping) and table S6).
