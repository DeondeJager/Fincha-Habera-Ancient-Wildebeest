# Admixture graph and f4-statistics in Admixtools2
For installation: https://uqrmaie1.github.io/admixtools/index.html
For tutorial: https://uqrmaie1.github.io/admixtools/articles/admixtools.html

## 1. Prepare input files
We used the `PACKEDANCESTRY` format in `Admixtools2`, which was produced from the `.bed`, `.bim`, `.fam` files made for the `EMU` [PCA](https://github.com/DeondeJager/Fincha-Habera-Ancient-Wildebeest/tree/main/Analyses/7-PCA) from the `n82_fossil` [pseudohaploid dataset](https://github.com/DeondeJager/Fincha-Habera-Ancient-Wildebeest/tree/main/Analyses/5-Pseudohaploid_genomes).

The `.bed`, `.bim`, `.fam` files were converted to the `PACKEDANCESTRY` format using the `convertf` function in `EIGENSOFT v8.0.0` (https://github.com/DReichLab/EIG).
 - This was done with the script `convertf_plink_to_admixtools2.sh`, which uses the file `convert_n82_fossil_BWD_HiC_autosomes_packedancestry.par`.

## 2. Run admixture graph analyses
Find the best-fititng graphs for each number of gene flow events from 0-6, with `find_graphs` implemented in the `admixtools2_qpGraph_subspecies_*.R` scripts, submitted as slurm cluster jobs by the `n82_fossil_BWD_HiC_admixtools2_subspecies_*.slurm` scripts.
 - These produce the R object files `opt_results_*.Rdata`, which contain a list of the best-fitting graph in each of the 100 `find_graphs` replicates implemented for each level of complexity (number of admixture events/admixture levels).

## 3. Process find_graphs output
The script `admixtools2_qpGraph_BWD_HiC_signif.R` reads the `winner_*.Rdata` files produced at step 2 and:
  
  i. Deduplicates the 100 graphs at each admxiture level based on their hash values,
  
  ii. Ranks the models from best to worst by their score (within admixture levels) and writes the results to `csv` files,
  
  iii. Plots the ranked models within admixture levels and saves the plots in `pdf` files,
  
  iv. Computes out-of-sample scores for comparing models within admixture levels and writes the results to `tsv` files,
  
  v. Calculates worst residuals as another measure by which to evaluate model fit and writes the results to `csv` files,
  
  vi. Performs bootstrap-resampled graph fits to statistically test whether different models within each winner group are significantly different or not and plots the p-values as a square matrix heatmap in a `pdf` file.

The script is submitted as a slurm cluster job with `n82_fossil_BWD_HiC_admixtools2_subspecies_signif.slurm`.

Because the last step takes very long it was submitted as a 14-day job on the cluster, which still did not finish and so had to be run for each admix level separately by changing "6" to "5" to "4", etc on lines 364-365 of `admixtools2_qpGraph_BWD_HiC_signif.R`.

For the final plots included in the manuscript, see the scripts `admixtools2_qpGraph_BWD_HiC_plots_final_Fig3_FigS14.R` and `admixtools2_qpGraph_BWD_HiC_plot_pMatrix_FigS14.R` in the [Figures](https://github.com/DeondeJager/Fincha-Habera-Ancient-Wildebeest/tree/main/Figures) folder.

## 4. f4-statistics
The f4-statistics analysis is in the Figures folder in the R script `Figures/Fig3_FigS14-Admixture_graphs/admixtools2_qpGraph_BWD_HiC_f4.R` and produces Fig. 3C.

Note that the `f2_results` in that folder contains the f2 statistics needed for the f4 analysis, pre-computed from the `packedancestrymap` file which is too large to share on GitHub.
