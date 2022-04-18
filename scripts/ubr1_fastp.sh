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
#SBATCH -o /home/albertf/mahlon/UBR1_RNA-seq_analysis/fastp/fastp_std_out_log
#SBATCH -e /home/albertf/mahlon/UBR1_RNA-seq_analysis/fastp/fastp_std_err_log

## -----
## setup

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

## [3] directory for fastp output
output_dir=${project_dir}fastp/output/
if [ ! -d ${output_dir} ]
then
    mkdir -v ${output_dir}
fi

## [4] sample listing that we use to read fastq files
## input_files=${project_dir}"test.txt"
input_files=${project_dir}"UBR1_sample_indices.txt"

## [5] logs
touch ${project_dir}fastp/fastp_std_{out,err}_log
out_log=${project_dir}fastp/fastp_std_out_log
err_log=${project_dir}fastp/fastp_std_err_log
echo "" > ${out_log}
echo "" > ${err_log}

## [6] fastp
fastp=~/.local/bin/fastp

## -----
## run fastp 
cd ${project_dir}fastp

while IFS=$'\t' read -r -a lineArray
do
echo "running fastp on sample: ${lineArray[0]}" >> ${out_log}
echo "running fastp genotype: ${lineArray[1]}" >> ${out_log}
${fastp} --in1 ${input_dir}${lineArray[1]}*R1*.fastq.gz \
--in2 ${input_dir}${lineArray[1]}*R2*.fastq.gz \
--out1 ${output_dir}${lineArray[1]}_fastp_out_R1.fastq.gz \
--out2 ${output_dir}${lineArray[1]}_fastp_out_R2.fastq.gz \
--cut_front --cut_tail --cut_mean_quality 15 \
--length_required 36 \
--dedup --detect_adapter_for_pe \
--trim_poly_g --trim_poly_x \
--report_title "${lineArray[1]}_fastp_report.html" \
--html "${lineArray[1]}_fastp_report.html"
done < ${input_files}
