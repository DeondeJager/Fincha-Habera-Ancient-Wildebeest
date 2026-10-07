# Convert pseudohaploid files from angsd to plink for use in emu

## Create plink tped and tfam files:
# get the converter from angsd website https://github.com/ANGSD/angsd/tree/master/misc/haplotoPlink.cpp
g++ -lz haploToPlink.cpp -o haplotoPlink.exe
# .exe file is also available with our data
chmod +x haplotoPlink.exe

## Navigate to directory with pseudohaploid dataset and convert
cd /path/to/dataset
~/programs/angsd/angsd-0.921/misc/haplotoPlink.exe n82_fossil_BWD_HiC_autosomes.haplo n82_fossil_BWD_HiC_autosomes

## Run Plink to make bed bim fam files used by emu
module load plink/1.9.0
plink --tfile n82_fossil_BWD_HiC_autosomes --make-bed --allow-extra-chr --autosome-num 28 --double-id --missing-genotype N --update-ids n82_fossil_BWD_HiC_recoded.txt --out n82_fossil_BWD_HiC_autosomes
