#!/bin/bash

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules and set up enviroment variables
module purge
module load getorganelle/1.7.7.0
READS=/projects/lorenzen/people/KUID/wildebeest/clean_data/Liu_2024
SEED=/projects/lorenzen/people/KUID/wildebeest/refgenomes
OUTDIR=/projects/lorenzen/people/KUID/wildebeest/GetOrganelle

# Navigate to working directory
cd $OUTDIR

## Without seed
#echo "Assembling SRR27955299 mitogenome WITHOUT seed..."
#get_organelle_from_reads.py -1 $READS/SRR27955299/SRR27955299_1.clean.fq.gz -2 $READS/SRR27955299/SRR27955299_2.clean.fq.gz --expected-max-size 16600 --expected-min-size 16350 -F animal_mt -o $OUTDIR/SRR27955299 -R 10 -t 4 --continue
#echo "Done assembling SRR27955299 mitogenome WITHOUT seed."
#date

## Initial assembly failed with "Disentangling timeout" and increasing --disentangle-time-limit did not help,
## thus, I used the get_organelle_from_assembly.py script to assemble the mitogenome from the assembly graph as follows:
#cd /projects/lorenzen/people/KUID/wildebeest/GetOrganelle/SRR27955299 # navigate to getorganelle folder for sample
#mkdir assemblyFromGraph # Make output directory
### Assemble mitogenome from K115 assembly graph with minmum depth settig to ignore contigs with low depth (likely contamination)
#get_organelle_from_assembly.py -F animal_mt -g extended_spades/K115/assembly_graph.fastg.extend-animal_mt.fastg -o assemblyFromGraph/ --min-depth 50 -t 4 --continue

## While the get_organelle_from_assembly.py worked, try again with settings that worked for Liu_59.sh after initial failure
## i.e. smaller word size (-w 100) and more rounds (-R 15)
get_organelle_from_reads.py -1 $READS/SRR27955299/SRR27955299_1.clean.fq.gz -2 $READS/SRR27955299/SRR27955299_2.clean.fq.gz -F animal_mt -o $OUTDIR/SRR27955299 -R 15 -t 4 -w 100 --overwrite

