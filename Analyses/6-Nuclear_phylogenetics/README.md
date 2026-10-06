# Nuclear phylogenetic analyses
Species tree reconstruction with the pseudohaploid datasets generated in `5-Pseudohaploid_genomes`:
- `*n82_fossil*`: Contains Fossil936 and the modern samples (excludes 277_damaged, but includes its undamaged bam, 277) + hartebeest outgroups.
- This applies to both the wildebeest BWD_HiC and goat reference genomes, as indicated in the file names where relevant.

>[!NOTE]
>While all files are provided in one folder here, each dataset should be run in its own directory, as some output files with the same name are produced and will thus will be overwritten.

## 1. Estimate phylogenies in 50 kb windows with phyml
The scripts `phyml_windows_pseudoHap_*.slurm` run Simon Martin’s `phyml_sliding_windows.py` [script](https://github.com/simonhmartin/genomics_general) to contruct the maximum likelihood phylogenies in 50kb windows from the pseudohaploid datasets generated in `5-Pseudohaploid_genomes`.
  - For the `phyml` tree estimation, they require an individuals file (`--indFile`) which lists all the names of the individual samples to include and must match the names in the pseudohaploid file. Here provided as `n82_fossil_all.txt`, `n81_modern_all.txt`, etc. 

## 2. Get species tree with ASTRAL
Thereafter, ASTRAL is used to estimate the species tree.
- This requires a mapping file (e.g. `n82_fossil_mapping_all.txt`), which lists the taxon groups and assigns the individuals to those groups, so ASTRAL knows which nodes to use for the quartet support analysis.
- See the ASTRAL tutorial [here](https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md) for more information.

## 3. Calculate quartet support values
ASTRAL may provide quartet scores already at the end of step 2, but there may also be no normalized quartet score provided with a message like "because of the existense of polytomies and to save time" in the `.log` output file.
 - I suspect I got this message because of many "polytomies" within subspecies and populations that the analysis was unable to resolve and not because the actual subspecies or species splits are polytomies.

In this case, see the `ASTRAL_score_trees.sh` file for how to proceed with scoring the trees to get both quartet scores and posterior probabilities for the nodes.
 - This requires the program `newick_utils`, which is a super useful program for working with phylogenies on the command-line: https://github.com/tjunier/newick_utils.
 - Newick utils is used to rename the tips from uninformative `ind0`-like names to the actual sample names using a renaming map file, like `n82_fossil_nw_rename_map.txt`.

## 4. Visualise trees
The final scored trees are provided here.

**BWD_HiC dataset**:
 - `n82_fossil.wOutgroup.50kb.ASTRAL.spTree.scored.rename.tre`: The final BWD_HiC tree with quartet scores.
 - `n82_fossil.wOutgroup.50kb.ASTRAL.spTree.scored.localPP.rename.tre`: The final BWD_HiC tree with with local posterior probabilities of the nodes.

**Goat dataset**:
 - `n82_fossil.Goat.wOutgroup.50kb.ASTRAL.spTree.scored.rename.tre`: The final goat tree with quartet scores.
 - `n82_fossil.Goat.wOutgroup.50kb.ASTRAL.spTree.scored.localPP.rename.tre`: The final BWD_HiC tree with with local posterior probabilities of the nodes.

The trees were visualised in [FigTree v1.4.4](http://tree.bio.ed.ac.uk/software/figtree/) and further edited for visual clarity in [Inkscape v1.2.1](https://inkscape.org).