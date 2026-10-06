#!/bin/bash

# Convert plink .bed, .bim, .fam (previously created for EMU PCA) to PACKEDANCESTRYMAP format for ADMIXTOOLS2

# Load eigensoft
module load eigensoft/8.0.0 # Contains the convertf program we need

# Do conversion
cd /path/to/plink/files
convertf -p convert_n82_fossil_BWD_HiC_autosomes_packedancestry.par