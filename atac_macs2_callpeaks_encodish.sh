#!/bin/bash
    
#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/macs_callPeaks_encodish.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/macs_callPeaks_encodish.err
#$ -q free64,som,asom,pub64 
#$ -pe make 16
#$ -R y 
#$ -t 1-6 
#$ -m beas            
#$ -ckpt blcr 

set -euxo pipefail

module load enthought_python/7.3.2
module load macs2/2.0.10
module load samtools/1.9

INPUT_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/analysis_encode/alignments/data_dedup_sorted
INPUT_LIST=${INPUT_DIR}/bam_list.txt

PEAK_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/analysis_encode/peak
RUNLOG_DIR=${PEAK_DIR}/runlogs

mkdir -p ${PEAK_DIR}
mkdir -p ${RUNLOG_DIR}

# for FILE in *.bam; do echo `pwd`/${FILE} >> sorted_nodup_bam_list.txt; done
INPUT=`head -n $SGE_TASK_ID $INPUT_LIST | tail -n 1`
INPUT_PREFIX=`head -n $SGE_TASK_ID $INPUT_LIST | tail -n 1 | cut -d. -f1 | cut -d/ -f 8 | cut -d_ -f 1`

# Capturing stout for each file
RUNLOG=${RUNLOG_DIR}/runlog_macs_peak_encodish_${SGE_TASK_ID}.txt
echo "Run by `whoami` on `date`" >> ${RUNLOG}
echo ${INPUT} >> ${RUNLOG}
echo ${INPUT_PREFIX} >> ${RUNLOG}

macs2 callpeak \
--treatment ${INPUT} \
--name ${PEAK_DIR}/${INPUT_PREFIX} \
--format BAM \
--nomodel \
--shift -100 \
--extsize 200 2>> ${RUNLOG}