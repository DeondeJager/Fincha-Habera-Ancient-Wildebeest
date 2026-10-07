# Get genome location of altitude target genes to extract them from bam files
module load samtools/1.23.1
cd /projects/lorenzen/people/KUID/wildebeest/refgenomes

## 1. Download BWD_HiC annotation files:
wget https://dnazoo.s3.wasabisys.com/Connochaetes_taurinus/BWD_HiC.fasta_v2.functional.gff3.gz
wget https://dnazoo.s3.wasabisys.com/Connochaetes_taurinus/BWD_HiC.fasta_v2.functional.proteins.fasta.gz
wget https://dnazoo.s3.wasabisys.com/Connochaetes_taurinus/BWD_HiC.fasta_v2.functional.transcripts.fasta.gz

## 2. Get genome location of EPAS1, as an example:
zgrep 'EPAS1' BWD_HiC.fasta_v2.functional.gff3.gz 
### The first result should be the one giving the location of the full gene - there will be several annotations with hits saying "match_part" to "HUMAN" or "MOUSE" - what you want is the top result, e.g.:
HiC_scaffold_9  maker   gene    28061790        28104362        .       +       .       ID=Connochaetes_taurinus011295;Name=Connochaetes_011295;Alias=maker-HiC_scaffold_9-exonerate_protein2genome-gene-281.1;Note=Similar to EPAS1: Endothelial PAS domain-containing protein 1 (Homo sapiens);

### Extract the gene region from the BWD_HiC reference genome
samtools faidx BWD_HiC.chr_Human_fossilMito.fasta HiC_scaffold_9:28061790-28104362 > BWD_HiC_EPAS1.fasta

### But when aligning this to the goat EPAS1 gene with CDS annotations, it was missing the first exon, which is ~52,000 np upstream.
### Thus, we extend the region 60 kb upstream to include this first exon:
samtools faidx BWD_HiC.chr_Human_fossilMito.fasta HiC_scaffold_9:28001790-28104362 > BWD_HiC_EPAS1.fasta
sed -i '1s/.*/>BWD_HiC_EPAS1/' BWD_HiC_EPAS1.fasta # Replace whole first line with new fasta header

## 3. Subset bam files to EPAS1 region
cd /projects/lorenzen/people/KUID/wildebeest/altitude_genes/EPAS1
### Fossil936
samtools view -b -o Fossil936_BWD_HiC_EPAS1.bam ../../paleomix/competitive_mapping/BWD_HiC/mapDamage/Fossil936.BWD_HiC-Human.merged.rescaled.rm.noX.bam HiC_scaffold_9:28001790-28104362
### Modern samples
for file in ../../bams/*BWD_HiC-Human.MD.MQ30.rm.noX.bam; do samtools view -b -o ${file/%.bam/_EPAS1.bam} $file HiC_scaffold_9:28001790-28104362; done

## Call fasta consensus sequences from the subset bams
### Fossil936
samtools consensus -o Fossil936_BWD_HiC_EPAS1.fasta --format fasta --mode bayesian --min-MQ 30 --min-BQ 20 --ambig --min-depth 1 Fossil936_BWD_HiC_EPAS1.bam
sed -i 's/HiC_scaffold_9/Fossil936_BWD_HiC_EPAS1/g' Fossil936_BWD_HiC_EPAS1.fasta # Rename fasta header
### Modern samples
for file in ../../bams/*BWD_HiC-Human.MD.MQ30.rm.noX.bam;  do samtools consensus -o ${file/%.bam/.fasta} --format fasta --mode bayesian --min-MQ 30 --min-BQ 20 --ambig --min-depth 1 $file; done
for file in *noX_EPAS1.fasta; do sed -i 's/-Human.MD.MQ30.rm.noX//' $file; done # Clean up fasta header by removing unnecessary text from filename string

