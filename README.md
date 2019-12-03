# epicPCR_analysis
Instructions and examples for how to analyze epicPCR data

# Prepare bioconda environment for epicPCR analysis

screen -S epicPCR
sinteractive
cd $WRKDIR/DONOTREMOVE/
module load bioconda/3
conda create -c bioconda --name epicPCR multiQC vsearch=2.6.0 fastx_toolkit cutadapt=1.10 pear=0.9.6 mothur=1.40.5 fastqc

## Before starting analysis activate the conda environment
source activate epicPCR


## Answer y when prompted

## Before starting analysis activate the conda environment
source activate epicPCR

## You can deactivate the environment like this
source deactivate epicPCR


# Adding environmental variables

### In addition to the bioinformatic programs we need to define some environmental variables which will make the analysis with the master script automatic. We need to tell the script where your SILVA database of SSUs is located and also what minimum and maximum length you want to use for trimming your reads. It is good to have a rather broad size range if you have several genes in your dataset. We also need to define the wanted threshold for clustering OTUs. 97% and 99% the are most commonly used.


## If you don't have a mothur version of the Silva database, download it from mothur website


mkdir SILVA132
cd SILVA132
wget https://mothur.org/w/images/3/32/Silva.nr_v132.tgz

tar zxvf Silva.nr_v132.tgz


## Set an environmental variable pointing to the Silva database. You need to change the path to your Silva database folder.

export SILVA_TAX="/wrk/YOURUSERNAME/DONOTREMOVE/PATHTOYOURSILVADBFOLDER/YOURSILVA.TAX"
export SILVA_ALN="/wrk/YOURUSERNAME/DONOTREMOVE/PATHTOYOURSILVADBFOLDER/YOURSILVA.ALN"


## Check that this works! If not you need to set the path again or change the name of your silva align and taxonomy files.
echo $SILVA_TAX
echo $SILVA_ALN



## Set variables for the minimum and maximum length of epicpCR products. This is dependent on YOUR OWN product lengths. Note that you might need to change these after you have done the analysis and looked at results. You can also filter out long and short sequences later so you don't need to set these too strict in the beginning.

export MIN_LEN="350"
export MAX_LEN="550"

## Set variables for the OTU clustering threshold

export OTU_TRH="0.99"

######################################################################

# Pre-analysis

### Before we start analysing we need to obtain some data to work with. We will use the data from Jenni Hultman's paper. Downloading files from ENA is a very useful tool for the future too so keep a hold of this script!

## Download files from Hultman et al. Epic project PRJEB23695 from ENA

##################################################################
## Open screen on taito
screen -S epic

## Log into interactive node
sinteractive

## Go to your DONOTREMOVE folder

cd $WRKDIR/DONOTREMOVE/

## Make directory myfiles
mkdir myfiles

## Go to your directory myfiles
cd myfiles

## Download the list of files to download
### This is obtained by going to the ENA website and searching with the project id that is available in the article. You can select which fields you want to get in the tabs by clicking. Here we want the submitted file names preserved so we have chosen the submitted_ftp format. 
wget -O list_of_files "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJEB23695&result=read_run&fields=submitted_ftp&download=txt"

## Modify list so that all samples are in one line
tr ';' '\n' < list_of_files > list_of_files2

mv -f list_of_files2 list_of_files

## Download files. This might take a while
while read file_name; do wget $file_name; done<list_of_files

## Uncompress. This might take a while
gunzip *gz 

## Check that you have all the files. Compare to the website source if needed.

## Good! Now you have all the raw data fastq files from Hultman et al. paper.


# Additional small files

### Next we need to make small files which are needed in the bioinformatic analysis. We need information about the 16S end primers, target gene linker (inside) and nested (outside/forward) primers and a list of simplified sample names we want to use further on. These are specific for YOUR OWN genes and YOUR OWN samples.


## Make file for 16S end primer
### The 785R needs to be in REVERSE COMPLEMENT orientation for this to work.
printf ">Illum_785R_1_RC\nGGATTAGATACCCNNGTAGTC\n>Illum_785R_2_RC\nGGATTAGATACCCNNGTAGTCt\n>Illum_785R_3_RC\nGGATTAGATACCCNNGTAGTCaga\n>Illum_785R_4_RC\nGGATTAGATACCCNNGTAGTCcactcag">785Rprimers.fasta



## Make file for primer mapping
### The file should not have a header and the first column has the gene name and the second column should have the nested primer or the outside/forward primer for epic functional gene and the third column the REVERSE COMPLEMENT of the linker or inside primer for epic functional gene.

## Jenni
printf "blaOXA\tTCGGTCTAAATGCGTGCCAT\tCTCATACTATGCTCAGCACA\ntetM\tGCAATTCTACTGATTTCTGC\tGCGGAAATTGTAATCAAACA\nqacEdelta\tGTATGCCGCATCTGAGGAAC\tTTTCTTTCTCTGGTTCTGAAATCCAT\nint1\tCGAAGTCGAGGCATTTCTGTC\tTCCTCGGTTTTCTGGAAGGC\n">map.txt

## Minna
printf "blaOXA\tTCGGTCTAAATGCGTGCCAT\tCTCATACTATGCTCAGCACA\ntetM\tGCAATTCTACTGATTTCTGC\tGCGGAAATTGTAATCAAACA\nstrB\tGTATGCCGCATCTGAGGAAC\tGTCAAACTGACTACGTCCAC\n">map.txt


## Make file for sample names
### The simplified sample name file should be in the same order as your sample reads are. You can check this manually after creating the sample_names file.
#Jenni
ls -r *R1_001.fastq | sed 's/A0..-//g' | awk -F '\\-Hultman' '{print $1}' | cut -c-5 | tr '[:lower:]' '[:upper:]' | tr '-' '_' > sample_names

## Minna
ls -r *R1_001.fastq | awk -F '-' '{print $1 "_" $2 "_" $3}' > sample_names


# The actual analysis of epicPCR data

### Now that we have the files we need with primer and sample name info and the raw data and our working environment set as we like to we can start the actual analysis.

### Remember to do this interactively in screen or as a batch job so you don't loose your work if your connection breaks.

## Join and trim files and remove 16S end primer

#!/bin/bash/

## Epic PCR analysis read pretreatment.

## Activate bioconda environment
module load bioconda/3
source activate epicPCR


## Analyze quality. This might take a while
mkdir -p fastqc
fastqc *fastq -o fastqc/.
multiqc fastqc/*.zip -n fastqc/multiqc_rawdata


## Print list of R1 reads and save the base name "read_base"
ls -tr *R1*fastq | sed 's/1_001.fastq//g' > read_base

## Remove 3' adapter from R1 and R2.
### (The reverse complement of the 16S (785R) primer's adapter is removed from R1 using option -a and the reverse complement of the ARG primer (F3) is removed from the R2) using option -A.

while read list; do cutadapt ./$list"1_001.fastq" ./$list"2_001.fastq" -a  AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A ATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT  -o ./$list"1_001_adapter_trimmed.fastq" -p $list"2_001_adapter_trimmed.fastq" ;done<read_base &> cutadapt_out_adapter_trimmed

fastqc *1_001_adapter_trimmed.fastq -o fastqc/.
fastqc *2_001_adapter_trimmed.fastq -o fastqc/.

multiqc fastqc/*_adapter_trimmed_fastqc.zip -n multiqc_adapter_trimmed

multiqc * -n multiqc_all_files

## Assemble reads with pear. This might take a while

## Jenni

while read list; do pear -y 150M -j 8 -f $list"1_001.fastq" -r $list"2_001.fastq" -o $list"12";done<read_base&>pear_out

## Minna

while read list; do pear -y 150M -j 8 -f $list"1_001_adapter_trimmed.fastq" -r $list"2_001_adapter_trimmed.fastq" -o $list"12";done<read_base&>pear_out

## List number of assembled files
ls -lt *.assembled.fastq | wc -l

## Remove not needed files
ls *discarded*

rm -f *discarded*

rm -f *unassembled*

## Analyze quality of merged reads. This might take a while

fastqc *R12.assembled.fastq -o fastqc/.

multiqc fastqc/*R12.assembled_fastqc.zip -n multiqc_R12.assembled

## Remove 16S end primer and do quality filtering.
### You can change the -q option based on the quality of the fastqc

while read name_base; do cutadapt ./$name_base*".assembled.fastq" -a file:785Rprimers.fasta -o $name_base"12_filtered.pair.fastq" --trimmed-only --max-n=5 -q 20 -m $MIN_LEN -M $MAX_LEN;done<read_base &> cutadapt_16S_out
less cutadapt_out

## Check with fastqc

fastqc *12_filtered.pair.fastq -o fastqc/. &>fastqc_out

ls -ltr *12_filtered.pair.fastq | wc -l

## Transform to fasta
while read name; do fastq_to_fasta -i $name"12_filtered.pair.fastq"  > $name"12_filtered.pair.fasta";done<read_base

## Make mapping file for renaming
paste sample_names read_base > name_mapping

## Rename the files
while read i
do
        arr=($i)
         mv -f ${arr[1]}12_filtered.pair.fasta ${arr[0]}_joined_assembled.fasta
done < name_mapping


## Rename for OTU table creation
while read i
do
        arr=($i)
sed "s/>@*/>barcodelabel=${arr[0]};read=/g"  ${arr[0]}_joined_assembled.fasta > ${arr[0]}_joined_assembled_renamed.fasta
done < sample_names

## Combine to one file
cat *_joined_assembled_renamed.fasta > all_joined_assembled.fasta

## Split the fasta file based on the genes using primer sequences

## Make folders for all target genes
while read i
do
        arr=($i)
	mkdir -p ${arr[0]}reads
done < map.txt

## Split scrip

### Extract reads starting with the forward primer (nested one)

while read i
 do
 arr=($i)

### Extract reads starting with the target gene forward primer (nested one) F3 (-g option)

cutadapt all_joined_assembled.fasta -g ${arr[1]} \
-o ${arr[0]}reads/filtered.${arr[0]}.fasta  --trimmed-only \
-O 10 -e 0.2 &> ${arr[0]}reads/filtered.${arr[0]}.log

### Add gene name to the sample name

sed "s/>barcodelabel=/>barcodelabel=${arr[0]}_/g" ${arr[0]}reads/filtered.${arr[0]}.fasta > ${arr[0]}reads/filtered.${arr[0]}.renamed.fasta


### From the reads which had the target gene primer, take the 16S region using the 16s primer 519F used in nested PCR (CAGCMGCCGCGGTAATWC) (-g option)
cutadapt ${arr[0]}reads/filtered.${arr[0]}.renamed.fasta \
-g CAGCMGCCGCGGTAATWC -o ${arr[0]}reads/filtered.${arr[0]}.16spart.fasta \
--trimmed-only -O 15  &>> ${arr[0]}reads/filtered.${arr[0]}.log

### Take the target gene region using the linker primer R1-F2' 's target gene region
cutadapt ${arr[0]}reads/filtered.${arr[0]}.renamed.fasta \
-a ${arr[2]} -o ${arr[0]}reads/filtered.${arr[0]}.${arr[0]}part.fasta \
 --trimmed-only -O 15 &>> ${arr[0]}reads/filtered.${arr[0]}.log  
 done < map.txt

#############################################################


## OTU clustering and making OTU tables using VSEARCH

### OTU clustering is optional. You can also only dereplicate the sequences and remove chimeras and classify all unique non-chimeric reads.

## Combine reads
cat *reads/filtered.*.16spart.fasta > all.filtered.16Sparts.fasta


## Dereplicate and discard singletons
vsearch --derep_fulllength all.filtered.16Sparts.fasta  --output all.unique.fasta --minuniquesize 2 --sizeout &> derep_out


## Make OTUs using vsearch's --cluster_fast
vsearch --cluster_fast all.unique.fasta --id $OTU_TRH --centroids 16S_OTUs.fasta --relabel OTU --uc 16Sclusters.uc &> clustering_out

## Map reads back to OTUs and produce OTU table
vsearch --usearch_global all.filtered.16Sparts.fasta  --db 16S_OTUs.fasta --strand plus --id $OTU_TRH --uc 16S_OTUtab.uc --otutabout 16S_OTUs.txt &> mapping_out

##########################################################################################################
##########################################################################################################
##########################################################################################################


## Classifying OTUs using mothur

mothur "#classify.seqs(fasta=16S_OTUs.fasta, reference=$SILVA_ALN,taxonomy=$SILVA_TAX, cutoff=60, processors=1, probs=F)"

## Modify taxonomy table file
sed 's/;/\t/gi' 16S_OTUs.nr_v132.wang.taxonomy > 16S.tax


## Optional
### Check the resistance gene against a database for example by blast

## Download resfinder and run blast

blastn -subject /wrk/parnanen/DONOTREMOVE/ARG_MGEdatabases/resfinder_FINAL.fa -query filtered.blaOXA.blaOXApart.fasta -outfmt 6 -out blaOXA_blast.out -max_target_seqs 1 

blastn -subject /wrk/parnanen/DONOTREMOVE/ARG_MGEdatabases/resfinder_FINAL.fa -query filtered.tetM.tetMpart.fasta -outfmt 6 -out tetM_blast.out -max_target_seqs 1

blastn -subject /wrk/parnanen/DONOTREMOVE/ARG_MGEdatabases/MGEs_FINAL.fa -query filtered.intI.intIpart.fasta -outfmt 6 -out intI_blast.out -max_target_seqs 1

blastn -subject /wrk/parnanen/DONOTREMOVE/ARG_MGEdatabases/MGEs_FINAL.fa -query filtered.qacEdelta.qacEdeltapart.fasta -outfmt 
