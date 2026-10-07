# Mapping aDNA data to goat reference
Map the Fossil936 and simulated SRR27955277_damaged reads to the goat reference genome.

This was not competitive mapping.

Instead, the reads mapping to the Fossil936 mitogenome, human nuclear genome, and human mitogenome in the BWD_HiC competitive mapping were filtered (removed) from the `reads.collapsed.gz` and `reads.collapsed.truncated.gz` fastq files produced by paleomix and the remaining reads were mapped to the goat nuclear genome only.
See below for an example of how this was done.

This essentially removes any potential human and NUMT contamination from the read pool and thus competitive mapping is not required.

- The `.slurm` files submit jobs that run paleomix with the settings in the `.yaml` files.

## Extract non-target reads from bam files
Non-target reads are those not mapping to the wildebeest BWD_HiC nuclear genome, i.e. mapping to human nuclear and mito genomes and to the Fossil936 mitogenome.

This is an example of the code. This was done for each Fossil936 library and the SRR27955277_damaged sample.

> **Note**: 
>
> This thus retains reads mapped to the wildebeest nuclear genome AND unmapped reads, so that we are not only keeping reads mapped to the wildebeest nuclear genome, as there might be endogenous reads in the unmapped pool that will map to the goat genome but not the wildebeest genome.

### Get list of non-target scaffolds
`-v` inverts the selection, and all wildebeest nuclear scaffolds start with "HiC_scaffold", so this selects all the scaffolds that do not start with "HiC_scaffold"
```
grep -v 'HiC_scaffold*' BWD_HiC.chr_Human_fossilMito.fasta.fai | cut -f1 > BWD_HiC.chr_Human_fossilMito.nonTarget.scaffolds.names.txt
```

### Extract reads mapping to the non-target scaffolds using samtools view and looping through the scaffold names list
The cut command extracts the read name (1st field), scaffold that it mapped to (2nd field), and read group (20th field) so that we can double-check it does not include any HiC_scaffold (wildebeest nuclear) scaffolds and we can see which library the read originated from (for removing the reads from those read pools later on).
```
module load samtools/1.17
rm SCR116.nonTarget.reads.txt # Need to remove any previously created output, as the command below doesn't overwrite but appends
while read -r line; do samtools view SCR116.BWD_HiC.chr_Human_fossilMito.bam $line | cut -f1,3,20 >> SCR116.nonTarget.reads.txt; done < BWD_HiC.chr_Human_fossilMito.nonTarget.scaffolds.names.txt
## Check scaffolds
cut -f2 SCR116.nonTarget.reads.txt | uniq
## Count reads
wc -l SCR116.nonTarget.reads.txt
```

### Filter reads
```
module load bbmap/39.80
## Navigate to the paleomix reads output folder
cd /projects/lorenzen/people/KUID/wildebeest/paleomix/competitive_mapping/allLibs/SCR116/SCR116/reads/SCR116/SCR116_DeepSeq/Lane_7/SCR116_EKDL230051128-1A_22FNFYLT3_L7_x.fq.gz/

## Filter reads using filterbyname.sh from bbmap
filterbyname.sh in=reads.collapsed.gz out=reads.collapsed.filter.gz names=SCR116.nonTarget.reads.txt include=f
```
These filtered reads were then mapped to the goat nuclear genome only with paleomix (see `.yaml` and `.slurm` files).
