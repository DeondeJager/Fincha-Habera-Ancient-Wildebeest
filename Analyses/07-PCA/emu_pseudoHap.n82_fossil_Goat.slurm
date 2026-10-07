#!/bin/bash
#SBATCH --job-name=emu_pseudoHap_Goat_n82_fossil    # Job name
#SBATCH --exclude=mjolnircomp05fl,mjolnircomp06fl # Exclude these two nodes, as they gave memory errors last time (they're older and have less memory (512GB) than the other nodes (1024GB).
#SBATCH --nodes=1                      # Run all processes on a single node
#SBATCH --ntasks=4                     # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G               # Memory required per CPU
#SBATCH --time=02:00:00                # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/emu/pseudohap/%x_%j.out    # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/emu/pseudohap/%x_%j.err     # Standard error log

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules/activate conda environments
## Initiate conda environment on the compute node first
source /home/KUID/.bashrc
### Activate package
conda deactivate
conda activate emu
module load plink/1.9.0
PSEUDOHAP=/projects/lorenzen/people/KUID/wildebeest/admixture_graph/Goat
REMOVE=/projects/lorenzen/people/KUID/wildebeest/plink/Goat
OUTDIR=/projects/lorenzen/people/KUID/wildebeest/emu/pseudohap/Goat
cd $OUTDIR

# Remove hartebeest from pseudohap plink files:
plink --allow-extra-chr --autosome-num 29 --double-id --remove $REMOVE/remove_Hartebeest.txt --bfile $PSEUDOHAP/n82_fossil_Goat_autosomes --out $OUTDIR/n82_fossil_Goat_autosomes_noHart --make-bed

# Run emu with pseudo-haploid file
## n82_fossil is with the fossil and the modern samples (excludes damaged SRR27955277, but includes its undamaged bam) + excludes hartebeest in this case for the pca
## Note: Transitions removed and maf=0.05 during creation of pseudohap file
emu --bfile n82_fossil_Goat_autosomes_noHart --eig 4 --threads 4 -o n82_fossil_Goat_autosomes_noHart.emu

# Remove hartebeest & black wildebeest from pseudohap plink files:
plink --allow-extra-chr --autosome-num 29 --double-id --remove $REMOVE/remove_Black.txt --bfile $PSEUDOHAP/n82_fossil_Goat_autosomes --out $OUTDIR/n72_fossil_Goat_autosomes_noBlack_noHart --make-bed

# Run emu with pseudo-haploid file
## n72_fossil is with the fossil and the modern blue samples (excludes damaged SRR27955277, but includes its undamaged bam, excludes black wildebeest) + excludes hartebeest in this case for the pca
## Note: Transitions removed and maf=0.05 during creation of pseudohap file
emu --bfile n72_fossil_Goat_autosomes_noBlack_noHart --eig 4 --threads 4 -o n72_fossil_Goat_autosomes_noBlack_noHart.emu

