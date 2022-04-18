#!/bin/bash
#PBS -l walltime=1:00:00,nodes=1:ppn=4,mem=32gb
#PBS -o /home/albertf/mahlon/illumina/2020.07.31_FPFA001_TDH3pr_dsRed_TFT_sorts/logs/
#PBS -e /home/albertf/mahlon/illumina/2020.07.31_FPFA001_TDH3pr_dsRed_TFT_sorts/error_logs/

#############
## USER INPUT
#############
## specify a folder for logs: 'PBS -o', above
## specify a folder for error logs: 'PBS -e', above
## trailing slashes at the end!!
#################
## END USER INPUT
#################

## load software; the '/xx' denotes a version 
bwaVersion='bwa/0.7.15'
module load ${bwaVersion}
samtoolsVersion='samtools/1.5'
module load ${samtoolsVersion}
module load bamtools

# the definitions above were for setting up. Now do this right, with getopt:
# see: 
# http://frodo.looijaard.name/project/getopt/misc/examples
# http://www.bahmanm.com/blogs/command-line-options-how-to-parse-in-bash-using-getopt

fileR1=${readDat}${fileRoot}${pairID1}${postfix}.${fileType}
fileR2=${readDat}${fileRoot}${pairID2}${postfix}.${fileType}
## find $fileR{1,2} 
## the loops below will throw an error if they don't exist

## 1. first loop -> ensure we have the first read file
if [ ! -e ${fileR1} ]
then
    echo "ERROR: File ${fileR1} not found. Exiting"; exit 1
fi

## 2. 2nd loop -> if we're in paired end mode, ensure 2nd read exists
## 3. if we pass both loops, align w/ samtools 
if [ "${PE}" = "PE" ]
then
    echo "Runing in PE mode"
    if [ ! -e ${fileR2} ]
    then
        echo "ERROR: File ${fileR2} not found in PE mode. Exiting"; exit 1
    fi
    bwa mem -t ${nThreads} ${genome} ${fileR1} ${fileR2} \
	| samtools sort -@${nThreads} -O BAM -o ${BWADat}"/"${fileRoot}sort.bam -
else
if [ "${PE}" = "SE" ]
then
    echo "Running in SE mode"
    bwa mem -t ${nThreads} ${genome} ${fileR1} \
	| samtools sort -@${nThreads} -O BAM -o ${BWADat}"/"${fileRoot}sort.bam -

else
    echo "ERROR: PE/SE mode not set correctly. Don't know what to do; Exiting"; exit 1
fi
fi

## the above will have either produced a sam or exited with a message

## get rid of mismatches
samtools view -q 30 ${BWADat}"/"${fileRoot}sort.bam \
    | grep 'XS:i:0' \
    | samtools view -b -T ${genome} - \
    > ${BWADat}"/"${fileRoot}sort_filtered.bam 

## get rid of PCR duplicates
samtools rmdup -S ${BWADat}"/"${fileRoot}sort_filtered.bam ${BWADat}"/"${fileRoot}sort_filtered_rmdup.bam 

## counting the coverage per SNPs
samtools mpileup -vu -t INFO/AD -l ${SNPs} -f ${genome} ${BWADat}"/"${fileRoot}sort_filtered_rmdup.bam \
	 > ${BWADat}"/"${fileRoot}sort_filtered_rmdup.vcf 

echo "done"
