# Process output fasta files (mitogenomes) of GetOrganelle to prep for annotation and rolling to same start position

# 1. Get summary statistics for all mitogenomes using summary_get_organelle_output.py
## Needs space-separated list of all directories that must be summarised
cd /projects/lorenzen/people/KUID/wildebeest/GetOrganelle
summary_get_organelle_output.py SRR27955210	SRR27955211 SRR27955212 SRR27955214 SRR27955215 SRR27955216 SRR27955217 SRR27955218 SRR27955219 SRR27955220 SRR27955221 SRR27955222 SRR27955223 SRR27955224 SRR27955225 SRR27955226 SRR27955227 SRR27955228 SRR27955229 SRR27955230 SRR27955231 SRR27955233 SRR27955234 SRR27955235 SRR27955236 SRR27955237 SRR27955238 SRR27955239 SRR27955240 SRR27955241 SRR27955242 SRR27955243 SRR27955244 SRR27955245 SRR27955246 SRR27955247 SRR27955248 SRR27955249 SRR27955250 SRR27955251 SRR27955252 SRR27955253 SRR27955254 SRR27955255 SRR27955256 SRR27955257 SRR27955258 SRR27955259 SRR27955260 SRR27955261 SRR27955262 SRR27955263 SRR27955264 SRR27955265 SRR27955266 SRR27955267 SRR27955268 SRR27955269 SRR27955270 SRR27955271 SRR27955272 SRR27955273 SRR27955274 SRR27955275 SRR27955276 SRR27955277 SRR27955278 SRR27955279 SRR27955280 SRR27955281 SRR27955283 SRR27955284 SRR27955285 SRR27955286 SRR27955287 SRR27955288 SRR27955289 SRR27955290 SRR27955291 SRR27955292 SRR27955294 SRR27955295 SRR27955296 SRR27955298 SRR27955299 SRR27955300 SRR27955301 SRR27955302 SRR27955303 SRR27955304 SRR27955305 SRR27955306 SRR27955307 SRR27955309 SRR27955310 SRR27955311 SRR27955312 SRR27955313 SRR27955314 SRR27955315 SRR27955318 SRR27955319 SRR27955320 SRR27955322 SRR27955324 SRR27955325 SRR27955326 SRR27955327 SRR27955328 SRR27955329 SRR27955330 SRR27955331 SRR27955332 SRR27955333 SRR27955335 SRR27955336 SRR27955337 SRR27955338 SRR27955339 SRR27955340 SRR27955341 SRR27955342 SRR27955344 SRR27955345 SRR27955346 SRR27955347 SRR27955348 SRR27955349 SRR27955350 SRR27955351 SRR27955352 SRR27955353 SRR27955354 -o GetOrganelle_Liu_stats.tsv

# 2. Copy and rename fasta files to match sample name (SRR)
cd /projects/lorenzen/people/KUID/wildebeest/GetOrganelle
## GetOrganelle output is in a subdirectory for each sample
while read -r line; do cp ${line}/*complete*.fasta ${line}/${line}.fasta; done < samples.txt
## samples.txt contains the name of each sample, one per line.

# 3. Rename fasta header to match the filename
## For inplace editing, you can use GNU awk (gawk), as follows:
## (Loops through the sample folders listed in samples.txt)
while read -r line; do
    cd "$line" || continue  # Change directory to the specified line (directory)
    gawk -i inplace -v INPLACE_SUFFIX=.bak '/^>/{print ">" substr(FILENAME, 1, length(FILENAME) - 6); next} 1' "${line}.fasta"
    cd ..  # Return to the parent directory
done < samples.txt
## Requires GNU awk >= v4.1.0, check with 'gawk --version' which version you have.
## The '-v INPLACE_SUFFIX=.bak' creates a backup file in case the inplace editing goes wrong.

# 4. Annotate using MITOS2
## This was just a test run. The actual annotation was done in a separate scripts here: /home/KUID/scripts/mitos2/Liu2024
module load mitos/2.1.9
cd SRR27955210
mkdir mitos2
runmitos.py -i SRR27955210.fasta --code 2 --outdir mitos2 --refdir /datasets/globe_databases/mitos/20201120/v0.3/ --refseqver refseq63m --debug --noplots
## This is successful (mostly - there are some plotting errors eventhough I suppressed plotting, but the main results are there)
## Extract the starting position of the trnF gene (first mtgene) to roll the mitogenome
STARTPOS=$(grep 'trnF' mitos2/result.gff | head -n1 | cut -f4) # Extract second column (start position) of trnF gene from bed file and save as environment variable
echo "${STARTPOS}" # check it worked
## Load mitobim to use the circules.py script to roll the mitogenome
module load mitobim/1.9.1
circules.py -f SRR27955210.fasta -n $STARTPOS
## The result is stored as a new fasta file called circular.rolled.${STARTPOS}.fasta
## Therefore, rename file and fasta header
SAMPLE=$(pwd | cut -f7 -d '/') # Get sample name
echo $SAMPLE
mv circular.rolled.${STARTPOS}.fasta ${SAMPLE}_rolled.fasta
## Rename fasta header to sample name (taken from filename)
#gawk -i inplace '/^>/{print ">" substr(FILENAME, 1, length(FILENAME) - 13); next} 1' ${SAMPLE}_rolled.fasta
## Or rename using the already stored sample name
awk -i inplace '/^>/{print ">" "'${SAMPLE}'"; next} 1' ${SAMPLE}_rolled.fasta
## Annotate the rolled mitogenome
mkdir mitos2_rolled
runmitos.py -i ${SAMPLE}_rolled.fasta --code 2 --outdir mitos2_rolled --refdir /datasets/globe_databases/mitos/20201120/v0.3/ --refseqver refseq63m --debug --noplots
## Copy gff output file to one with the same name as the sample
cp mitos2_rolled/result.gff ./${SAMPLE}_rolled.gff

