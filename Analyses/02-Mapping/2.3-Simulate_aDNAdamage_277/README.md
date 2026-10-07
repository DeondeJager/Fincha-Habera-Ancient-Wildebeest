# Simulate aDNA damage on sample SRR27955277 and map reads

## 1. Simulate aDNA damage
- `aDNA_mutate_modern.sh`: Shortens reads with `fastp` and adds aDNA with `tapas`. Author: Michael V. Westbury
- Uses the `Mutationprofile.tab` file for simulating aDNA damage. 
  - This file was produced from the Fossil936 mapDamage `misincorporation.txt` output file from the BWD-HiC bam file.
- See `tapas` docs for details: https://mlell.github.io/tapas/03_read-mutation.html

## 2. Map simulated aDNA reads

Paleomix was used to map the reads to the wildebeest BWD-HiC competitive genome (see 2.1-Fossil936_mapping) and goat reference genome.
- `SRR27955277_paleomix_compMap_BWD_HiC.yaml`: Contains paleomix parameters; called by slurm script `Wildebeest_SRR27955277_compMap_BWD_HiC.slurm`
- `SRR27955277_paleomix_goat.yaml`: Contains paleomix parameters; called by slurm script `Wildebeest_SRR27955277_goat.slurm`
