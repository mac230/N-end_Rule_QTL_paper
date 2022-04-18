#!/bin/bash -l
#SBATCH --time=35:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH -p ram256g
#SBATCH -n 2
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahlon@umn.edu
#SBATCH -o /home/albertf/mahlon/UBR1_RNA-seq_analysis/fastp/kallisto_std_out_log
#SBATCH -e /home/albertf/mahlon/UBR1_RNA-seq_analysis/fastp/kallisto_std_err_log

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
echo ${project_dir}

## [2] directory w/ fastq files
input_dir=${project_dir}fastp/output/
ls -alh ${input_dir}

## [3] directory for kallisto output
output_dir=${project_dir}fastp/kallisto_out/
if [ ! -d ${output_dir} ]
then
    mkdir -v ${output_dir}
fi
ls -alh ${output_dir}


## [4] sample listing that we use to read fastq files
## input_files=${project_dir}"test.txt"
input_files=${project_dir}"UBR1_sample_indices.txt"
head ${input_files}


## [5] kallisto command - using version 46.1, which is 
## not available on MSI
kallisto=~/kallisto/kallisto
${kallisto} version

## [6] log files
touch ${project_dir}/fastp/kallisto_std_{out,err}_log
out_log=${project_dir}/fastp/kallisto_std_out_log


## -----
## build the index
## '-i' specifies the index to be built
${kallisto} index -i "cerevisiae_transcripts.idx" \
${project_dir}/kallisto/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa


## -----
## run kallisto
## IFS = "internal field separator"; how bash
## recognizes word boundaries when splitting
## sequences of character strings; here we set
## it to tab to match the table we supply
## see: https://unix.stackexchange.com/questions/18886/why-is-while-ifs-read-used-so-often-instead-of-ifs-while-read
## for more information.  Here, the input is
## specified at the end of the loop
## ls -l ${input_dir}

## kallisto parameters
## quant: run the quantification algorithm
## -i: specifies the index (get from the kallisto page)
## -t: threads, here 8
## -l: est'd avg. fragment length
## -s: est'd std dev of fragment length, here 30 bp
## -b: n. bootstrap samples, here 100
## -o: where to write output to 
## NOTE: avg. frag. length and length sd automatically estimated
## for paired-end data, so can drop these parameters
cd ${project_dir}/fastp/

while IFS=$'\t' read -r -a lineArray
do
echo "$( date +%Y.%m.%d_%T ) running kallisto on sample: ${lineArray[1]}" > ${out_log}
${kallisto} quant \
--index=s_cer_transcriptome.idx \
-t 8 --rf-stranded -b 100  \
-o ${output_dir}${lineArray[1]}"_kallisto_out" \
${input_dir}${lineArray[1]}*R1*.fastq.gz \
${input_dir}${lineArray[1]}*R2*.fastq.gz \
--genomebam --gtf Saccharomyces_cerevisiae.R64-1-1.96.gtf \
--chromosomes chrLengths_ensemblFormat.txt
done < ${input_files}
