#!/bin/bash
    
#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/encode_coverage.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/encode_coverage.err
#$ -q free64,som,asom,pub64 
#$ -pe make 16
#$ -R y 
#$ -t 1-6 
#$ -m as            
#$ -ckpt blcr

set -euxo pipefail

module load bedtools/2.25.0

# Input files
BED_LIST=/som/sborrego/201810_ATACSEQ_MB468_R8/analysis_encode/peak/narrowPeak_list.txt
BAM_LIST=/som/sborrego/201810_ATACSEQ_MB468_R8/analysis_encode/alignments/data_dedup_sorted/bam_list.txt

# Taskarray file names
BED_INPUT=`head -n $SGE_TASK_ID $BED_LIST | tail -n 1`
BAM_INPUT=`head -n $SGE_TASK_ID $BAM_LIST | tail -n 1`
COVER_PREFIX=`head -n $SGE_TASK_ID $BED_LIST | tail -n 1 | cut -d. -f1 | cut -d_ -f5 | cut -d/ -f3`

# Output directory
COVER_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/analysis_encode/peak/coverage
COVER_OUTPUT=${COVER_DIR}/${COVER_PREFIX}_coverage.txt

mkdir -p $COVER_DIR

RUNLOG_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/analysis_encode/peak/runlogs
RUNLOG=${RUNLOG_DIR}/runlog_coverage_${SGE_TASK_ID}.txt
echo "Run by `whoami` on `date`" >> ${RUNLOG}
echo ${BAM_INPUT} >> ${RUNLOG}
echo ${BED_INPUT} >> ${RUNLOG}
echo ${COVER_OUTPUT} >> ${RUNLOG}

bedtools coverage \
-abam ${BAM_INPUT} \
-b ${BED_INPUT} > ${COVER_OUTPUT} 2>> ${RUNLOG}