#!/bin/bash

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load module
module load fastp/0.23.2

# Navigate to working directory
cd /projects/lorenzen/people/KUID/wildebeest/raw_data/Liu_2024

# Insert commands here:
/home/KUID/scripts/readsQC/fastp_wrapper_modern_Liu2024_wildebeest.sh \
 -d /projects/lorenzen/people/KUID/wildebeest \
 -r raw_data/Liu_2024 -s SRR27955234 -t 10

