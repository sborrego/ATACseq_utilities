#!/bin/bash
    
#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/encode_tagAlign.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/encode_tagAlign.err
#$ -q free64,som,asom,pub64 
#$ -pe make 16
#$ -R y 
#$ -t 1-6 
#$ -m beas            
#$ -ckpt blcr

set -euxo pipefail

module load bedtools/2.22
module load picard-tools/1.96

# Paths to Picard programs
SORTSAM=/data/apps/picard-tools/1.96/SortSam.jar

# Input files
BAM_LIST=/som/sborrego/201810_ATACSEQ_MB468_R8/analysis_encode/alignments/data_dedup_sorted/bam_list

# Output directories
TAG_ALIGN_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/analysis_encode/alignments/tagAlign
NAMESORT_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/analysis_encode/alignments/data_dedup_name_sorted

mkdir -p ${NAMESORT_DIR}
mkdir -p ${TAG_ALIGN_DIR}

# Taskarray selection of files
BAM_INPUT=`head -n $SGE_TASK_ID $BAM_LIST | tail -n 1`
BAM_OUTPUT=`head -n $SGE_TASK_ID $BAM_LIST | tail -n 1 | cut -d. -f1 | cut -d/ -f 8 | cut -d_ -f 1`

# Creating virtual SE file containing both read pairs
FINAL_TAG_ALIGN_FILE=${TAG_ALIGN_DIR}/${BAM_OUTPUT}

RUNLOG_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/analysis_encode/alignments/runlogs
RUNLOG=${RUNLOG_DIR}/runlog_tagAlign_{SGE_TASK_ID}.txt
echo "Run by `whoami` on `date`" >> ${RUNLOG}
echo ${BAM_INPUT} >> ${RUNLOG}
echo ${BAM_OUTPUT} >> ${RUNLOG}

# Create tagAlign File
bedtools bamtobed \
-i ${BAM_INPUT} \
awk 'BEGIN{OFS="\t"}{$4="N"; $5="1000"; print $0}' \
| gzip -c > ${FINAL_TAG_ALIGN_FILE}.PE2SE.tagAlign.gz 2>> ${RUNLOG}

# Name sorting deduped files for BEDPE file
java -Xmx2g -jar ${SORTSAM} \
	INPUT=${BAM_INPUT} \
	OUTPUT=${NAMESORT_DIR}/${BAM_OUTPUT}.dedup.namesorted.bam \
	SORT_ORDER=queryname \
	VALIDATION_STRINGENCY=LENIENT 2>> ${RUNLOG}

# Create BEDPE File

bedtools bamtobed \
-bedpe \
-mate1 \
-i ${NAMESORT_DIR}/${BAM_OUTPUT}.dedup.namesorted.bam \
| gzip -c > ${FINAL_TAG_ALIGN_FILE}.bedpe.gz 2>> ${RUNLOG}

# Subsample tagAlign file
# Restrict to one read end per pair for CC analysis

NREADS=15000000
SUBSAMPLED_TA_FILE="${FINAL_TAG_ALIGN_FILE}.filt.nodup.sample.$((NREADS /
1000000)).MATE1.tagAlign.gz"

zcat ${FINAL_TAG_ALIGN_FILE}.bedpe.gz \
| grep -v “chrM” \
| shuf -n ${NREADS} \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"N","1000",$9}' \
| gzip -c > ${SUBSAMPLED_TA_FILE} 2>> ${RUNLOG}
