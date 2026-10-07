# Modern wildebeest and hartebeest outgroup read mapping

## Process raw reads
Remove any adapters and quality trim the modern reads downloaded from the SRA.

The script `fastp_wrapper_modern_Liu2024_wildebeest.sh` contains the actual `fastp` parameters implemented.

The folder `Liu202_fastp` contains the per-sample bash scripts that call the above `fastp` script.

The `fastp_wrapper_Liu2024_wildebeest.batch.slurm` slurm script submits these as an array job on the cluster.

## Map reads to BWD_HiC-Human and goat reference genomes
The modern samples were mapped against the same "competitive" wildebeest-human-nuclear-mito reference genome as Fossil936, to make downstream processing and analyses consistent and easier.

See the subfolder `Liu2024_mapping` for the per-sample bash scripts, which are submitted as an array jobs by the `Liu2024_wildebeest.BWD_HiC.batch.slurm` and `Liu2024_wildebeest.Goat.batch.slurm` scripts.

Those scripts call the script `MapWGS_nuclear_Liu2024_wildebeest.sh`, which is where the actual `bwa`, `samtools`, `picard` etc. commands can be found.
