#!/bin/bash
    
#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/openRegions.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/openRegions.err
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
TAG_DIR=${PEAK_CALLING_DIR}/tag_directories
PEAKS_DIR=${PEAK_CALLING_DIR}/peaks
PEAK_200=${PEAKS_DIR}/peak_200
PEAK_500=${PEAKS_DIR}/peak_500
PEAK_MERGE=${PEAKS_DIR}/peak_merge
PEAK_FINAL=${PEAKS_DIR}/peak_final
RUNLOG_DIR=${PEAK_CALLING_DIR}/runlogs

mkdir -p ${TAG_DIR}
mkdir -p ${PEAKS_DIR}
mkdir -p ${PEAK_200}
mkdir -p ${PEAK_500}
mkdir -p ${PEAK_MERGE}
mkdir -p ${RUNLOG_DIR}

# for FILE in *.bam; do echo `pwd`/${FILE} >> sorted_nodup_bam_list.txt; done
INPUT=`head -n $SGE_TASK_ID $INPUT_LIST | tail -n 1`
INPUT_PREFIX=`head -n $SGE_TASK_ID $INPUT_LIST | tail -n 1 | cut -d/ -f7 | cut -d. -f1`

# Capturing stout for each file
RUNLOG=${RUNLOG_DIR}/runlog_peakCalling_${SGE_TASK_ID}.txt
echo "Run by `whoami` on `date`" >> ${RUNLOG}
echo ${INPUT} >> ${RUNLOG}
echo ${INPUT_PREFIX} >> ${RUNLOG}

# Make one tag directory per bam file
mkdir -p ${TAG_DIR}/${INPUT_PREFIX}
makeTagDirectory ${TAG_DIR}/${INPUT_PREFIX} ${INPUT} 2>> ${RUNLOG}

# Find open regions aroud 200 
mkdir -p ${PEAK_200}/${INPUT_PREFIX}

findPeaks ${TAG_DIR}/${INPUT_PREFIX} \
	-o ${PEAK_200}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak200.txt \
	-style factor \
	-size 200 \
	-fdr 0.01 2>> ${RUNLOG}

# Find open regions aroud 500 
mkdir -p ${PEAK_500}/${INPUT_PREFIX}

findPeaks ${TAG_DIR}/${INPUT_PREFIX} \
	-o ${PEAK_500}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak500.txt \
	-style factor \
	-size 500 \
	-fdr 0.01 2>> ${RUNLOG}

# Merge files with open regions of 200bp and 500bp
bedops -merge \
${PEAK_200}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak200.txt \
${PEAK_500}/${INPUT_PREFIX}/${INPUT_PREFIX}_peak500.txt > ${PEAK_MERGE}/${INPUT_PREFIX}_bedops.merge.txt 2>> ${RUNLOG}

# Sort merged files
sortBed \
-i ${PEAK_MERGE}/${INPUT_PREFIX}_bedops.merge.txt > ${PEAK_MERGE}/${INPUT_PREFIX}_bedops.merge.sorted.txt 2>> ${RUNLOG}

# Merge single sorted file - results in final peak list
bedtools merge \
-i ${PEAK_MERGE}/${INPUT_PREFIX}_bedops.merge.sorted.txt > ${PEAK_FINAL}/${INPUT_PREFIX}_finalpeaks.txt 2>> ${RUNLOG}
