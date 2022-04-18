#!/bin/sh

## once this finishes, make sure all files succeeded
## e.g. by counting the number of vcfs in the alignments folder:
## ls -alh ${BWADat}*vcf | wc -l
## compare this to the number of fastq files in the data dir: 
## ls -alh ${readDat}*R1*.gz | wc -l

#############
## USER INPUT
#############
## trailing "/" at the end of directories!!
## your data location here
readDat='/home/albertf/data_release/umgc/nextseq/200730_NB551498_0052_AH5MVTAFX2/Albert_Project_017/'

## folder where alignments will go
BWADat='/home/albertf/mahlon/illumina/2020.07.31_FPFA001_TDH3pr_dsRed_TFT_sorts/alignments/'

## specify where the logs go in 'mapFilterCount_4Cluster_FPFA001.sh' (lines 4,5)
## 'mapscript' = the script we'll run in the 'for' loop below
mapScript='/home/albertf/mahlon/illumina/2020.07.31_FPFA001_TDH3pr_dsRed_TFT_sorts/scripts/mapFilterCount_4Cluster_FPFA001.sh'
#################
## END USER INPUT
#################

## these usually won't change, but double check:
genome='/home/albertf/shared/genomes/sacCer3.fa'
## be sure the SNP set is two columns for the reference genome
## and that it has the correct line breaks (safe in xcode once, if necessary)
SNPs='/home/albertf/shared/SNPSets/SNPs_Maggie_170809_BY_positions.txt'
nThreads=24
fileType='fastq.gz'
PE='PE'
pairID1='R1'
pairID2='R2'
postfix='_001'

## need to be careful that the underscores in front of ${pairID1} match
fileRoots=( $(ls ${readDat} | grep "${pairID1}" | sed "s/${pairID1}.*//") )

## make sure this worked properly, e.g.: 
## printf '%s\n' ${fileRoots[*]} | wc -l
## ls -alh ${readDat}*R1*.gz | wc -l

for thisRoot in ${fileRoots[4]}; do

    qsub -v fileRoot=${thisRoot},PE=${PE},nThreads=${nThreads},readDat=${readDat},BWADat=${BWADat},genome=${genome},SNPs=${SNPs},pairID1=${pairID1},pairID2=${pairID2},postfix=${postfix},fileType=${fileType} ${mapScript}

## test by hand for a single sample by replacing '${thisRoot}' w/ '${fileRoots[4]}', e.g.
## run as a single line starting at 'qsub -v', above (i.e., don't start the 'for' loop)

done
