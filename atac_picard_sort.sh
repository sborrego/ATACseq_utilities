#!/bin/bash
    
#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/samsort.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/samsort.err
#$ -q free64,som,asom,pub64 
#$ -pe make 16
#$ -R y 
#$ -t 1-6 
#$ -m beas            
#$ -ckpt blcr  

set -euxo pipefail

module load java/1.7
module load picard-tools/1.96
module load samtools/1.0

SORTSAM=/data/apps/picard-tools/1.96/SortSam.jar
BAM_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/alignments/data_trim_clip_noSORT/chr_only_noSORT_alignments
BAM_LIST=${BAM_DIR}/bam_list.txt

BAM_INPUT=${BAM_DIR}/`head -n $SGE_TASK_ID $BAM_LIST | tail -n 1`
BAM_OUTPUT=`head -n $SGE_TASK_ID $BAM_LIST | tail -n 1 | cut -d_ -f1`

SORTED_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/alignments/data_chr_trim_clip_SORTED

mkdir -p ${SORTED_DIR}

RUNLOG=${SORTED_DIR}/runlog_samsort_${SGE_TASK_ID}.txt
echo "Run by `whoami` on `date`" > ${RUNLOG}
echo ${BAM_INPUT} > ${RUNLOG}
echo ${BAM_OUTPUT} > ${RUNLOG}

java -Xmx2g -jar SORTSAM \
INPUT=${BAM_INPUT} \
OUTPUT=${SORTED_DIR}/${BAM_OUTPUT}.sorted.bam \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=LENIENT 2>> ${RUNLOG}