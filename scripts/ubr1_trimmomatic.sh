#!/bin/bash -l
#SBATCH --time=15:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH -p ram256g
#SBATCH -n 2
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahlon@umn.edu
#SBATCH -o /home/albertf/mahlon/UBR1_RNA-seq_analysis/trimmomatic/std_out_log
#SBATCH -e /home/albertf/mahlon/UBR1_RNA-seq_analysis/trimmomatic/std_err_log

## -----
## setup

## variables
## [1] project directory
## determine which system we're on and set accordingly
if [ "$(hostname)" = mahlon-linux ]
then
    project_dir=~/emacs/N-end_Rule_QTL_paper/UBR1_RNA-seq_analysis/
else
    project_dir=~/UBR1_RNA-seq_analysis/
fi

## [2] directory w/ fastq files
input_dir=/home/albertf/data_release/umgc/nextseq/210728_VH00601_2_AAAGV5CHV/Albert_Project_023/

## [3] directory for trimmomatic output
output_dir=${project_dir}trimmomatic/
if [ ! -d ${output_dir} ]
then
    mkdir -v ${output_dir}
fi

## [4] sample listing that we use to read fastq files
## input_files=${project_dir}"test.txt"
input_files=${project_dir}"UBR1_sample_indices.txt"


## -----
## load software
module load trimmomatic/0.39


## -----
## run trimmomatic
## IFS = "internal field separator"; how bash
## recognizes word boundaries when splitting
## sequences of character strings; here we set
## it to tab to match the table we supply
## see: https://unix.stackexchange.com/questions/18886/why-is-while-ifs-read-used-so-often-instead-of-ifs-while-read
## for more information.  Here, the input is
## specified at the end of the loop
# ls -l ${input_dir}
TRIMMOMATIC=/panfs/roc/msisoft/trimmomatic/0.39/

while IFS=$'\t' read -r -a lineArray
do
echo "sample: ${lineArray[0]}"
echo "genotype: ${lineArray[1]}"
java -jar ${TRIMMOMATIC}trimmomatic-0.39.jar PE -threads 14 \
${input_dir}${lineArray[1]}*R1*.fastq.gz \
${input_dir}${lineArray[1]}*R2*.fastq.gz \
-baseout "${output_dir}${lineArray[1]}.fastq.gz" \
ILLUMINACLIP:${TRIMMOMATIC}adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done < $input_files

## parameters:
## LEADING: remove bases from the beginning w/ quality < 3
## TRAILING: remove bases from the end w/ quality < 3
## SLIDINGWINDOW: cut once the average quality within a 4 bp win < 15
## MINLEN: remove reads w/ length < 36

## ls -l ${output_dir}

## seems there's multiple adapter files, may need to use
## the other one if this doesn't work.  see:
## https://www.biostars.org/p/250425/
## cat $TRIMMOMATIC/adapters/TruSeq3-PE-2.fa
## cat $TRIMMOMATIC/adapters/TruSeq3-PE.fa
