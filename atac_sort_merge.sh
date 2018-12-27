#!/bin/bash
    
#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/mergeNsort.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/mergeNsort.err
#$ -q free64,som,asom,pub64 
#$ -pe make 16
#$ -R y 
#$ -t 1-6 
#$ -m beas            
#$ -ckpt blcr  

set -euxo pipefail

module load samtools/1.9
module load bedops/2.4.14
module load bedtools/2.23.0
module load homer/4.7

INPUT_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/alignments/data_final_181213
INPUT_LIST=${INPUT_DIR}/sorted_nodup_bam_list.txt

PEAK_CALLING_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/peak_calling
PEAKS_DIR=${PEAK_CALLING_DIR}/peaks
PEAK_200=${PEAKS_DIR}/peak_200
PEAK_500=${PEAKS_DIR}/peak_500
PEAK_MERGE=${PEAKS_DIR}/peak_merge
PEAK_FINAL=${PEAKS_DIR}/peak_final
RUNLOG_DIR=${PEAK_CALLING_DIR}/runlogs

mkdir -p ${PEAK_200}
mkdir -p ${PEAK_500}
mkdir -p ${PEAK_MERGE}
mkdir -p ${PEAK_FINAL}
mkdir -p ${RUNLOG_DIR}

# for FILE in *.bam; do echo `pwd`/${FILE} >> sorted_nodup_bam_list.txt; done
INPUT=`head -n $SGE_TASK_ID $INPUT_LIST | tail -n 1`
INPUT_PREFIX=`head -n $SGE_TASK_ID $INPUT_LIST | tail -n 1 | cut -d/ -f7 | cut -d. -f1`

# Capturing stout for each file
RUNLOG=${RUNLOG_DIR}/runlog_sortNmerge_${SGE_TASK_ID}.txt
echo "Run by `whoami` on `date`" >> ${RUNLOG}
echo ${INPUT_PREFIX} >> ${RUNLOG}

# Removing header from peak files
cat ${PEAK_200}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak200.txt | sed '/^\#/d' > ${PEAK_200}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak200_headless.txt
cat ${PEAK_500}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak500.txt | sed '/^\#/d' > ${PEAK_500}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak500_headless.txt

# Convert peak files to bed format
pos2bed.pl ${PEAK_200}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak200_headless.txt > ${PEAK_200}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak200.not_fixed.bed 2>> ${RUNLOG}
pos2bed.pl ${PEAK_500}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak500_headless.txt > ${PEAK_500}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak500.not_fixed.bed 2>> ${RUNLOG}

# Fixing column numbers that may be reversed
awk '{OFS="\t"} {if ($3<$2) print $1,$3,$2 ; else print $0}' ${PEAK_200}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak200.not_fixed.bed > ${PEAK_200}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak200.bed
awk '{OFS="\t"} {if ($3<$2) print $1,$3,$2 ; else print $0}' ${PEAK_500}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak500.not_fixed.bed > ${PEAK_500}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak500.bed

# Sort and merge bed files
sort -k1,1 -k2,2n \
${PEAK_200}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak200.bed \
${PEAK_500}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak500.bed \
| mergeBed -c 4,5 -o sum,collapse -i - > ${PEAK_MERGE}/${INPUT_PREFIX}.merge.bed 2>> ${RUNLOG}

# Merge single sorted file - results in final peak list
bedtools merge -c 4,5 -o sum,collapse \
-i ${PEAK_MERGE}/${INPUT_PREFIX}.merge.bed > ${PEAK_FINAL}/${INPUT_PREFIX}.finalpeaks.bed 2>> ${RUNLOG}

