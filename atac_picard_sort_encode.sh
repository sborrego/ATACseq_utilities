#!/bin/bash
    
#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/encode_dedup_samsort.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/encode_dedup_samsort.err
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

# Paths to Picard programs
MARK_DUPS=/data/apps/picard-tools/1.96/MarkDuplicates.jar
SORTSAM=/data/apps/picard-tools/1.96/SortSam.jar

# Input files
BAM_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/analysis_encode/alignments/sorted_alignments
BAM_LIST=${BAM_DIR}/bam_list.txt

# Output directories
DEDUP_SORTED_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/alignments/data_dedup_sorted
DEDUP_METRICS_DIR=${DEDUP_SORTED_DIR}/metrics

mkdir -p ${DEDUP_SORTED_DIR}
mkdir -p ${DEDUP_METRICS_DIR}

# Building file names
BAM_INPUT=`head -n $SGE_TASK_ID $BAM_LIST | tail -n 1`
BAM_OUTPUT=`head -n $SGE_TASK_ID $BAM_LIST | tail -n 1 | cut -d/ -f8 | cut -d. -f1`
BAM_FINAL=${DEDUP_SORTED_DIR}/${BAM_OUTPUT}

RUNLOG_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/analysis_encode/alignments/runlogs
RUNLOG=${RUNLOG_DIR}/runlog_dedup_sort_${SGE_TASK_ID}.txt
echo "Run by `whoami` on `date`" >> ${RUNLOG}
echo ${BAM_INPUT} >> ${RUNLOG}
echo ${BAM_OUTPUT} >> ${RUNLOG}

FLAG=${RUNLOG_DIR}/runlog_dedup_sort_${SGE_TASK_ID}.flagstat
echo "Run by `whoami` on `date`" >> ${FLAG}

java -Xmx2g -jar ${MARK_DUPS} \
	INPUT=${BAM_INPUT} \
	OUTPUT=${BAM_FINAL}.dedup.bam \
	METRICS_FILE=${DEDUP_METRICS_DIR}/${BAM_FINAL}.metrics.txt \
	REMOVE_DUPLICATES=true 2>> ${RUNLOG}

samtools flagstat ${BAM_FINAL}.dedup.bam >> ${FLAG}

java -Xmx2g -jar ${SORTSAM} \
	INPUT=${BAM_FINAL}.dedup.bam \
	OUTPUT=${DEDUP_SORTED_DIR}/${BAM_FINAL}.dedup.sort.bam \
	SORT_ORDER=coordinate \
	VALIDATION_STRINGENCY=LENIENT 2>> ${RUNLOG}

samtools index ${DEDUP_SORTED_DIR}/${BAM_FINAL}.dedup.sorted.bam
samtools flagstat ${DEDUP_SORTED_DIR}/${BAM_FINAL}.dedup.sorted.bam >> ${FLAG}