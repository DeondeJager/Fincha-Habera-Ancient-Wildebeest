# Estimate heterozygosity for Fossil936 and contemporary samples
Based on the results of the heterozygosity testing done in `4.1-Het_tests_277`, we estimated genome-wide heterozygosity for all samples using `ANGSD` with the `baq 1` option, as well as no transitions `-noTrans 1`, amongst other filters (see scripts).

## Fossil936
- `angsd_het_Fossil936_BWD_HiC.slurm`
  - This produces a `.ml` file with three colums and a row for each "window" of 10 million sites specified by `-nSites` in the `realSFS` command.

The proportion of heterozygous sites in each window is then calculated from the `.ml` output file by dividing the second column (heterozygous sites) by the sum of all three columns (where columns one and three represent the sites homozygous for the reference and alternate alleles, respectively).

## Contemporary samples
See folder `Liu2024_het` for per-sample commands.
The scripts are submitted as slurm array jobs by `angsd_het_Liu2024_BWD_HiC.batch.slurm`.

Of note `-setMinDepth` and `-setMaxDepth` were individually set for each sample based on the median and standard deviation of its coverage, as calculated by Picard’s CollectWgsMetrics in the mapping script ([2.2-Modern_samples_mapping/MapWGS_nuclear_Liu2024_wildebeest.sh](https://github.com/DeondeJager/Fincha-Habera-Ancient-Wildebeest/blob/main/Analyses/2-Mapping/2.2-Modern_samples_mapping/MapWGS_nuclear_Liu2024_wildebeest.sh), and see table S6 for these values).

## Plotting and statistics
The heterozygosity of the modern samples was calculated from the `*.ml` files for the modern samples in the same R script used for plotting the results (see [Figures](https://github.com/DeondeJager/Fincha-Habera-Ancient-Wildebeest/tree/main/Figures)).

The statistical comparisons of mean heterozygosity of Fossil936 and the modern subspecies are also conducted in the R plotting script, entitled `plot_het_BWD_HiC_Fig2F.R` in the Figures folder. 
