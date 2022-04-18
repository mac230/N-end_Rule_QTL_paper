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
#SBATCH -o /home/albertf/mahlon/UBR1_RNA-seq_analysis/fastp/rseqc/rseqc_std_err_log
#SBATCH -e /home/albertf/mahlon/UBR1_RNA-seq_analysis/fastp/rseqc/rseqc_std_out_log


## -----
## setup
## http://rseqc.sourceforge.net/#download-rseqc
## suggests installing using the python installer:
## which pip
## pip install RSeQC
## note that the scripts get installed in:
## ~/.local/bin/ after the above pip command

## have to load R so that geneBody_coverage.py
## can create the output pdf 
module load R/3.6.0

base_dir=/home/albertf/mahlon/UBR1_RNA-seq_analysis/fastp/
ls -alh ${base_dir}
mkdir -v ${base_dir}rseqc/
rseq_dir=${base_dir}rseqc/
## ls -alh ${rseq_dir}

## -----
## rsync as needed
## move files from home to msi
## rsync -avzn ~/emacs/N-end_Rule_QTL_paper/UBR1_RNA-seq_analysis/rseqc/
## mahlon@mesabi.msi.umn.edu:/home/albertf/mahlon/UBR1_RNA-seq_analysis/kallisto/rseqc

## move files from msi to home
## rsync -avzn mahlon@mesabi.msi.umn.edu:/home/albertf/mahlon/UBR1_RNA-seq_analysis/kallisto/rseqc
## ~/emacs/N-end_Rule_QTL_paper/UBR1_RNA-seq_analysis/rseqc/


## -----
## logs
r_std_out=${rseq_dir}rseqc_std_out_log
r_std_err=${rseq_dir}rseqc_std_err_log
touch ${r_std_out}
touch ${r_std_err}
echo "" > ${r_std_out}
echo "" > ${r_std_err}
cat ${r_std_out}
cat ${r_std_err}


## -----
## run rseqc scripts
## variables
gene_bod=~/.local/bin/geneBody_coverage.py
## ${gene_bod} --version
gene_ref=/home/albertf/mahlon/UBR1_RNA-seq_analysis/kallisto/rseqc/ensemblGenes_fromUCSC_mod.bed
## head ${gene_ref}

## -----
## 'tin.py'
## measures the 'transcript integrity number' (TIN):
## TIN, an algorithm that is able to measure RNA integrity at transcript level.
## calculates a score (0 <= TIN <= 100) for each expressed transcript
## rgeneseqc.sourceforge.net/#download-rseqc
tin=~/.local/bin/tin.py
${tin} --version

names='wild-type mutant'

for i in {1..5};
do
for n in ${names};
do
input=${n}_0${i}
echo "running tin.py on sample:" ${input} >> ${r_std_out}
## write the output to the appropriate dir
cd ${base_dir}kallisto_out/${input}_kallisto_out
echo "$( date ) running geneBody_coverage.py on sample:" ${input} >> ${r_std_out}
${tin} -i ${base_dir}kallisto_out/${input}_kallisto_out/${input}*pseudoalignments.bam -r ${gene_ref} \
>> ${r_std_out} 2>> ${r_std_err}
done
done


## -----
## 'geneBody_coverage.py':
## calculate the RNA-seq reads coverage over gene body.
## rgeneseqc.sourceforge.net/#download-rseqc
## the output here is a pdf of a line graph 
## showing the coverage across transcripts 

cd ${rseq_dir}

names='wild-type mutant'

for i in {1..5};
do
for n in ${names};
do
input=${n}_0${i}
echo "running geneBody_coverage.py on sample:" ${input} >> ${r_std_out}
${gene_bod} -i ${base_dir}/kallisto_out/${input}_kallisto_out/${input}*pseudoalignments.bam \
-r ${gene_ref} -o ${input} >> ${r_std_out} 2>> ${r_std_err}
done
done
