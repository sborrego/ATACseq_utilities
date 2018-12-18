#!/bin/bash
    
#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/markDuplicates.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/markDuplicates.err
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

# Paths to piccard programs
MARK_DUPS=/data/apps/picard-tools/1.96/MarkDuplicates.jar
INSERT_SIZE=/data/apps/picard-tools/1.96/CollectInsertSizeMetrics.jar

# Input files
BAM_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/alignments/data_chr_trim_clip_SORTED
BAM_LIST=${BAM_DIR}/sorted_bam_list.txt

# Output directories
NO_DUPS_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/alignments/data_final_181213
HIST_DIR=${NO_DUPS_DIR}/histograms

mkdir -p ${NO_DUPS_DIR}
mkdir -p ${HIST_DIR}

# Building file names
BAM_INPUT=${BAM_DIR}/`head -n $SGE_TASK_ID $BAM_LIST | tail -n 1`
BAM_OUTPUT=`head -n $SGE_TASK_ID $BAM_LIST | tail -n 1 | cut -d. -f1`
BAM_FINAL=${NO_DUPS_DIR}/${BAM_OUTPUT}.sorted.nodup.bam

RUNLOG=${NO_DUPS_DIR}/runlog_markDuplicates_${SGE_TASK_ID}.txt
echo "Run by `whoami` on `date`" >> ${RUNLOG}
echo ${BAM_INPUT} >> ${RUNLOG}
echo ${BAM_OUTPUT} >> ${RUNLOG}

FLAG=${NO_DUPS_DIR}/picard_markDuplicates_stats.flagstat
echo "Run by `whoami` on `date`" >> ${FLAG}

java -Xmx2g -jar ${MARK_DUPS} \
	INPUT=${BAM_INPUT} \
	OUTPUT=${BAM_FINAL} \
	METRICS_FILE=${NO_DUPS_DIR}/${BAM_OUTPUT}.marked_dup_metrics.txt \
	REMOVE_DUPLICATES=true 2>> ${RUNLOG}

samtools index ${BAM_FINAL} 
samtools flagstat ${BAM_FINAL} >> ${FLAG}

java -Xmx2g -jar ${INSERT_SIZE} \
      INPUT=${BAM_FINAL} \
      OUTPUT=${HIST_DIR}/${BAM_OUTPUT}.insert_size_metrics.txt \
      HISTOGRAM_FILE=${HIST_DIR}/${BAM_OUTPUT}.insert_size_histogram.pdf \
      MINIMUM_PCT=0.5 2>> ${RUNLOG}

