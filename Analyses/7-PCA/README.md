# Principal component analysis
PCA with EMU: https://github.com/Rosemeis/emu.

## 1. Convert pseudohaploid datasets to PLINK files
See the script `convert_pseudohap_toPlink.sh`, which:
 - Converts the pseudohaploid dataset to PLINK `.tped` and `.tfam` files using the `haploToPlink` subprogram from [ANGSD](https://github.com/ANGSD/angsd/tree/master/misc), which needed to be compiled first (instructions in script, but executable also provided here: `haplotoPlink.exe`).
 - Converts the `.tped` and `.tfam` files to `.bed`, `.bim`, `.fam` files need by `EMU` using `PLINK v1.9.0`.
   - This step also uses the file `n82_fossil_BWD_HiC_recoded.txt` to rename the individuals for easier plotting downstream.
   - These `*recoded.txt` files were used for both the BWD_HiC and goat datasets, as they contained the same individuals in the same order.

## 2. Perform PCA with EMU
### BWD_HiC dataset
- `emu_pseudoHap.n82_fossil_BWD_HiC.slurm`:
  - Includes analyses of `n82_fossil` and `n81_modern` datasets (see [5-Pseudohaploid_genomes](https://github.com/DeondeJager/Fincha-Habera-Ancient-Wildebeest/blob/main/Analyses/5-Pseudohaploid_genomes) for explanations).
  - Includes analyses of the above datasets with all wildebeest and with blue wildebeest only (i.e. black wildebeest removed).
  - Produces output files used to generate Fig. 2C and fig. S10 (see [Figures](https://github.com/DeondeJager/Fincha-Habera-Ancient-Wildebeest/tree/main/Figures) folder).

### Goat dataset
- `emu_pseudoHap.n82_fossil_Goat.slurm`:
  - Includes analyses of `n82_fossil` dataset with all wildebeest and with blue wildebeest only (i.e. black wildebeest removed).
  - Produces files used to generate fig. S8C (see [Figures](https://github.com/DeondeJager/Fincha-Habera-Ancient-Wildebeest/tree/main/Figures) folder).

Both of these scripts use the `remvoe_Hartbeest.txt` and `remove_Black.txt` files to remove these samples from the pseudohaploid datasets prior to PCA analysis, where appropriate.
