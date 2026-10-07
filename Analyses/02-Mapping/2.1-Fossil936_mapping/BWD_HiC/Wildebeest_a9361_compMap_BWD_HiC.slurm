#!/bin/bash
#SBATCH --job-name=a9361_paleomix_compMap_BWD_HiC         # Job name
#SBATCH --nodes=1                       # Run all processes on a single node
#SBATCH --ntasks=1                      # Run on a single task
#SBATCH --cpus-per-task=20               # Number of CPU cores per task
#SBATCH --mem-per-cpu=20G                # Memory required per CPU
#SBATCH --time=14-00:00:00                 # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/paleomix/competitive_mapping/BWD_HiC/a9361/a9361_paleomix_compMap_BWD_HiC_%j.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/paleomix/competitive_mapping/BWD_HiC/a9361/a9361_paleomix_compMap_BWD_HiC_%j.err         # Standard error log

# The %j in the --output & --error lines tells SLURM to substitute the job ID in the name of the output file
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts 
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script
hostname; date

echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules
module load paleomix/1.3.8

# Navigate to working directory where yaml file is
cd /projects/lorenzen/people/KUID/wildebeest/paleomix/competitive_mapping/BWD_HiC/a9361

# Execute paleomix
paleomix bam run a9361_paleomix_compMap_BWD_HiC.yaml --max-threads 20 --adapterremoval-max-threads 20 --bwa-max-threads 20

# Remove old summary file (otherwise it will be appended and not overwritten with the new data)
rm a9361_paleomix_compMap_BWD_HiC_summary.txt

# Combine results into summary table
awk 'BEGIN{second=0;one="";two="";three="";header="";line=""}{if($0!~/^#/ && $0!=""){if($0!~/^Target/){if($1!=one || $2!=two || $3!=three){if(header!="" && second==0){print "Target\tSample\tLibrary\t"header;second=1};if(line!=""){print one"\t"two"\t"three"\t"line;line=""}};if(second==0){split($0,h,"#");header=header h[2] "\t"}line=line $5 "\t";one=$1;two=$2;three=$3}}}END{print one"\t"two"\t"three"\t"line}' *summary | grep "\*\s*\*\s*" >> a9361_paleomix_compMap_BWD_HiC_summary.txt

