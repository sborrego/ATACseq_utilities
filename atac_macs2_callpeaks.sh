#!/bin/bash
    
#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/macs_callPeaks.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/macs_callPeaks.err
#$ -q free64,som,asom,pub64 
#$ -pe make 16
#$ -R y 
#$ -t 1-6 
#$ -m beas            
#$ -ckpt blcr 

set -euxo pipefail

module load enthought_python/7.3.2
module load macs2/2.0.10

INPUT_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/alignments/data_final_181213
INPUT_LIST=${INPUT_DIR}/sorted_nodup_bam_list.txt

PEAK_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/alignments/macs_peaks
RUNLOG_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/alignments/macs_peaks/runlogs

mkdir -p ${PEAK_DIR}
mkdir -p ${RUNLOG_DIR}

# for FILE in *.bam; do echo `pwd`/${FILE} >> sorted_nodup_bam_list.txt; done
INPUT=`head -n $SGE_TASK_ID $INPUT_LIST | tail -n 1`
INPUT_PREFIX=`head -n $SGE_TASK_ID $INPUT_LIST | tail -n 1 | cut -d/ -f7 | cut -d. -f1`

# Capturing stout for each file
RUNLOG=${RUNLOG_DIR}/runlog_macs_peakCalling_${SGE_TASK_ID}.txt
echo "Run by `whoami` on `date`" >> ${RUNLOG}
echo ${INPUT} >> ${RUNLOG}
echo ${INPUT_PREFIX} >> ${RUNLOG}


macs2 callpeak
--treatment ${INPUT} \
--name ${INPUT_PREFIX} \
--format AUTO \
--nomodel \
--shift -100 \
--extsize 200 2>> ${RUNLOG}