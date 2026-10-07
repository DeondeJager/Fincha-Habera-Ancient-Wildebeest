#!/bin/bash

# This is here to help calculate duration of analysis - see end of file for more
SECONDS=0

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "[Project: Fincha Habera wildebeest]"
   echo
   echo "This script is a wrapper around fastp for QC and trimming of paired-end Illumina reads from modern DNA. The script assumes that the fastp executable is available in one way or another (i.e. the correct conda environment needs to be activated, or the module loaded (if it exists), or it must be in your $PATH."
   echo
   echo "[USAGE]: fastp_wrapper_modern_Liu2024_wildebeest.sh [-d|h|r|t]"
   echo "Options:"
   echo "d: Path to parent data directory. This path must not contain any spaces or special characters. The raw data directory (option 'r' below) must be inside this directory. Default is /projects/lorenzen/people/KUID/wildebeest/. Within this parent folder must be 'raw_data' and 'clean_data' subfolders, but 'clean_data' will be created automatically if it does not exist."
   echo "h: Print this Help."
   echo "r: Raw sequencing data parent folder name (e.g. raw_data). Must not contain spaces or special characters. In the 'raw_data' subfolder must be another subfolder for each sample, which then contains the forward and reverse reads in *1.fq.gz *2.fq.gz format. It is assumed there is one library per modern sample, but there can be multiple lanes per library. If there were more libraries per sample, then each should have their own subfolder with the *.fq.gz files - uncomment 'MULTIPLE LIBS' sections of the code when this is the case and comment-out the 'SINGLE LIB' sections."
   echo "s: Sample name. There should be a folder for each sample within the raw_data folder (-r) that contains the reads for that sample."
   echo "t: Number of threads to use for fastp. This should be equal to or less than the number of CPUs requested via SLURM. Default is 1."
   echo
}

############################################################
############################################################
# Main program                                             #
############################################################
############################################################

# Set variables
Data=/projects/lorenzen/people/KUID/palaeobovids
Threads=1

############################################################
# Process the input options. Add options as needed.        #
############################################################

# Get the options
# Note: if an option is followed by a colon ":" in the "while getops..." line, 
# then it requires an argument, but if not (as with "h") then no argument is required
while getopts d:hr:s:t: option; do
   case "$option" in
      d) # Enter path to parent data directory
         Data=$OPTARG;;
      h) # display Help
         Help
         exit;;
      r) # Enter name of raw data folder/directory
         Rawdata=$OPTARG;;
      s) # Enter name of sample
         Sample=$OPTARG;;
      t) # Enter the number of threads to use for fastp
         Threads=$OPTARG;;
      \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done

# For sanity checks and testing
echo
echo "Sanity check of directory and settings:"
echo
echo "The parent data folder is (option d): $Data"
echo
echo "The raw data is in this here (option r): $Data/$Rawdata"
echo
echo "The sample being processed is (option s): $Sample"
echo
echo "The number of threads being used is (option t): $Threads"
echo
echo "The reads must be in here: $Data/$Rawdata/$Sample"
echo
#sleep 5m

# Actual commands
# Navigate to data directory and
# create a list of all the read files/libraries to be processed. i.e. Make a list of the files (subdirectories), one per line, in the raw_data folder and save to file. This file will then be used in the fastp command to loop through each readgroup/library.
cd $Data/$Rawdata/$Sample

#--- SINGLE LIB START ---#
ls -1 *_1.fq.gz > reads.txt # List forward reads files, one per line. Each *_1.fq.gz must have a corresponding *_2.fq.gz file
#--- SINGLE LIB END ---#

#--- MULTIPLE LIBS START ---#
#ls -l | grep ^d  | sed 's:.*\ ::g' > dir.list # Uncomment this line
#--- MULTIPLE LIBS END ---#

# Make output directory/ies for fastp if they do not already exist. If they already exist, then mkdir will give an error telling you this, which can be ignored.
mkdir $Data/clean_data
mkdir $Data/clean_data/$Sample

#--- MULTIPLE LIBS START ---#
#for line in $(cat $"dir.list"); do
#    mkdir $Data/$Rawdata/clean_data/$line
#done
#--- MULTIPLE LIBS END ---#

#--- SINGLE LIB START ---#
# Run fastp for each readgroup
for line in $(cat $"reads.txt"); do
 fastp --in1 $Data/$Rawdata/$Sample/${line} \
 --in2 $Data/$Rawdata/$Sample/${line/_1/_2} \
 --out1 $Data/clean_data/Liu_2024/$Sample/${line%_1.fq.gz}_1.clean.fq.gz \
 --out2 $Data/clean_data/Liu_2024/$Sample/${line%_1.fq.gz}_2.clean.fq.gz \
 --unpaired1 $Data/clean_data/Liu_2024/$Sample/${line%_1.fq.gz}_1.unpaired.fq.gz \
 --unpaired2 $Data/clean_data/Liu_2024/$Sample/${line%_1.fq.gz}_2.unpaired.fq.gz \
 --failed_out $Data/clean_data/Liu_2024/$Sample/${line%_1.fq.gz}.failed_out.fq.gz \
 --compression 4 \
 --correction \
 --overlap_len_require 30 --overlap_diff_limit 5 --overlap_diff_percent_limit 20 \
 --qualified_quality_phred 20 --unqualified_percent_limit 40 \
 --average_qual 20 \
 --length_required 30 \
 --trim_poly_g \
 --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
 --overrepresentation_analysis --overrepresentation_sampling 100 \
 --html $Data/clean_data/Liu_2024/$Sample/fastp.${line%_1.fq.gz}.html \
 --json $Data/clean_data/Liu_2024/$Sample/fastp.${line%_1.fq.gz}.json \
 --report_title "${line%_1.fq.gz} fastp report" \
 --thread $Threads
done
#--- SINGLE LIB END ---#

#--- MULTIPLE LIBS START ---#
# Run fastp for each set of reads (i.e. sequencing lane). Uncomment all lines below
#for line in $(cat $"dir.list"); do
# fastp --in1 $Data/$Rawdata/$Sample/$line/${line}*1.fq.gz \
# --in2 $Data/$Rawdata/$Sample/$line/${line}*2.fq.gz \
# --out1 $Data/clean_data/$Sample/${line}_1.clean.fq.gz \
# --out2 $Data/clean_data/$Sample/${line}_2.clean.fq.gz \
# --unpaired1 $Data/clean_data/$Sample/$line.unpaired1.fq.gz \
# --unpaired2 $Data/clean_data/$Sample/$line.unpaired2.fq.gz \
# --failed_out $Data/clean_data/$Sample/$line.failed_out.fq.gz \
# --compression 4 \ 
# --correction \
# --overlap_len_require 30 --overlap_diff_limit 5 --overlap_diff_percent_limit 20 \
# --qualified_quality_phred 20 --unqualified_percent_limit 40 \
# --average_qual 20 \
# --length_required 30 \
# --trim_poly_g \
# --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
# --overrepresentation_analysis --overrepresentation_sampling 100 \
# --html $Data/clean_data/$Sample/fastp.$line.html \
# --report_title "${line} fastp report" \
# --thread $Threads
#done
#--- MULTIPLE LIBS END ---#

#wait

# How long did the analysis take?
echo
duration=$SECONDS
echo "The analysis took $(($duration / 3600)) hours $(($duration / 60)) minutes and $(($duration % 60)) seconds to complete."
