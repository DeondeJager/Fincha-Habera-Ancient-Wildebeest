# Genetic sexing of Fossil936

## 1. Determine which HiC_scaffold represents the X chromosome

### 1.1. Find which scaffolds have the amelogenin gene annotated
The amelogenin (*AMEL*) gene is a sex-linked gene in mammals, found on the X and Y chromosomes.
- Download the file `BWD_HiC.fasta_v2.functional.gff3.gz` from the DNA Zoo [assembly page](https://dnazoo.s3.wasabisys.com/index.html?prefix=Connochaetes_taurinus/) of BWD-HiC
- Extract lines containing the word "AMEL":
```
zgrep 'AMEL' BWD_HiC.fasta_v2.functional.gff3.gz > BWD_HiC.fasta_v2.functional.AMEL.gff3
```
- This file (provided here) shows annotations for *AMELX* and *AMELY* on both HiC_scaffold_3 and HiC_scaffold_22.

### 1.2. Perform synteny analysis of HiC_scaffolds with bovine X and Y chromosomes
For this analysis I used the `SatsumaSynteny.sh` script from the [SeXY](https://github.com/andreidae/SeXY) pipeline, also included here.

I used NC_037357.1 (Cattle_X) and CM001061.2 (Cattle_Y) as the sex chromosome target sequences, concatenated into a single fasta file `RefX_RefY.fasta`, 
and the BWD-HiC assembly `BWD_HiC.chr.fasta` as the query sequence, where I retained only the chromosome-level HiC_scaffolds 1-29.

- Run the SatsumaSynteny analysis with `SatsumaSynteny.slurm`. 
- The `satsuma_summary.chained.out` file (included here) contains the syntenic regions identified between the target and query sequences.
- In this file, we see syntenic regions between 19 of the HiC_scaffolds and either the cattle chromosome X or chromosome Y. 
However, we note that (i) most of these matches are short, (ii) HiC_scaffold_22 is absent (identified as potential sex chromosome in step 1.1), and (iii) HiC_scaffold_3 has the most and longest syntenic regions to both the cattle chromosome X and Y sequences.

Thus, as HiC_scaffold_3 is the only one identified in both 1.1 and 1.2, it most likely represents the wildebeest X chromosome.

To confirm this, I conducted one more analysis, where I compared the mean depth of coverage of HiC_scaffold_3 and the rest of the HiC_scaffolds of a known male and female wildebeest from the modern samples in step 1.3.

### 1.3. Confirm HiC_scaffold_3 as X chromosome with modern samples
- I mapped reads from a modern wildebeest known male (SRR27955227) and known female (SRR27955224) [Liu et al. 2024](https://doi.org/10.1038/s41467-024-47015-y) to the BWD_HiC reference genome using the same script provided in the folder `2-Mapping`.
- Then calculated the mean depth of coverage 10 times at 1,000,000 randomly sampled sites across all HiC_scaffolds using the `CoverageCalculation.sh` from the SeXY pipeline.
- This produces the files `SRR27955224_female_HiC_scaffolds_coverage.txt` and `SRR27955227_female_HiC_scaffolds_coverage.txt`, which report the 10 replicates of coverage calculation for each HiC_scaffold.
  - These files were then used to produce fig. S5 (see [Figures](https://github.com/DeondeJager/Fincha-Habera-Ancient-Wildebeest/tree/main/Figures) folder).

These results show that HiC_scaffold_3 has approximately half the coverage of the rest of the scaffolds in the male sample, whereas the female sample has equal coverage across all chromosomes, **confirming HiC_scaffold_3 as the wildebeest X chromosome**.

## 2. Determine genetic sex of Fossil936
I could then determine the sex of Fossil936 using the relative coverage between HiC_scaffold_3 and the rest of the scaffolds (autosomes).

### 2.1. Re-map subset of reads
Because HiC_scaffold_3 had been previously been masked as the X chromosome in the bam file, 
I had to re-map Fossil936 reads to calculate the relative coverage of this scaffold and the rest (autosome).

For this, I used only a subset of the reads , to save time and resources, and SeXY is well-suited to low coverage data, so this is not a problem.
- Reads used: a9361_HiSeq (ERR16741334) and a9362_HiSeq (ERR16741336)
- Paleomix settings in `Fossil936_HiSeq_paleomix_compMap_BWD_HiC_SeXY.yaml`, submitted as a slurm job by `Wildebeest_HiSeq_compMap_BWD_HiC_SeXY.slurm`

### 2.2. Calculate coverage of sex chromosome and autosomes
`Fossil936_CalculateCoverage.slurm` calls `CoverageCalculation.sh` from the SeXY pipeline, which calculates the mean depth of coverage 10 times at 1,000,000 randomly sampled sites.
 - Also uses the files `HiC_scaffold_3.bed` and `HiC_autosomes.bed`.
 - Produces the file `Fossil936.BWD_HiC-Human_ratios.txt`, where the 10 chrX:Autosomes ratios can be found. 
 - This was plotted in R to produce Fig. 1C (see the [Figures](https://github.com/DeondeJager/Fincha-Habera-Ancient-Wildebeest/tree/main/Figures) folder).
