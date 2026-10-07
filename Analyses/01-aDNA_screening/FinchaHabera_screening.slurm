#!/bin/bash
#SBATCH --job-name=FinchaHabera_screening_mito_panel.slurm         # Job name
#SBATCH --exclude=mjolnircomp05fl,mjolnircomp06fl # Exclude these two nodes, as they're a bit older and slower than the others.	
#SBATCH --nodes=1                       # Run all processes on a single node
#SBATCH --ntasks=1                      # Run on a single task
#SBATCH --cpus-per-task=10               # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G                # Memory required per CPU
#SBATCH --time=24:00:00                 # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/KUID/wildebeest/paleomix/%x_%j.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/KUID/wildebeest/paleomix/%x_%j.err         # Standard error log

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
conda activate paleomix # This is version 1.3.9

# Navigate to working directory where yaml file is
cd /projects/lorenzen/people/KUID/wildebeest/paleomix/screening

# Execute paleomix
echo "$(date) Mapping Fincha Habera libraries to mito panel..."
paleomix bam run paleomix_20201126_Pool1_mito_panel --max-threads 10 --adapterremoval-max-threads 10 --bwa-max-threads 10
echo "$(date) Done mapping Fincha Habera libraries to mito panel."

# Create summary
## Remove previous file if it exists, otherwise it will be appended
rm FinchaHabera_mito_panel_summary.txt
## Aggregate results
awk 'BEGIN{second=0;one="";two="";three="";header="";line=""}{if($0!~/^#/ && $0!=""){if($0!~/^Target/){if($1!=one || $2!=two || $3!=three){if(header!="" && second==0){print "Target\tSample\tLibrary\t"header;second=1};if(line!=""){print one"\t"two"\t"three"\t"line;line=""}};if(second==0){split($0,h,"#");header=header h[2] "\t"}line=line $5 "\t";one=$1;two=$2;three=$3}}}END{print one"\t"two"\t"three"\t"line}' *summary | grep "\*\s*\*\s*" >> FinchaHabera_mito_panel_summary.txt

