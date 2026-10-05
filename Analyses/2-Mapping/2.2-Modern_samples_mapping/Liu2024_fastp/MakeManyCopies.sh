# Make copies of original script and replace the "_1" with the next number
# Written with help from Google's Gemini Advanced
for i in {2..133}; do cp Liu_1.sh Liu_$i.sh; done

# Replace the SRR number from Liu_1.sh with a new number as read from SRR accession list provided
i=2 # Set counter to 2, so that Liu_1.sh is ignored
# Then replace the string - make sure that you have deleted "SRR27955210" from the input file
while read -r new_str; do sed -i "s/SRR27955210/$new_str/g" Liu_$i.sh; i=$((i+1)); done < SRR_Acc_List_Included_132.txt

