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

MARK_DUPS=/data/apps/picard-tools/1.96/MarkDuplicates.jar

BAM_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/alignments/data_chr_trim_clip_SORTED
BAM_LIST=${BAM_DIR}/sorted_bam_list.txt
NO_DUPS_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/alignments/data_final_181213

mkdir -p ${NO_DUPS_DIR}

BAM_INPUT=${BAM_DIR}/`head -n $SGE_TASK_ID $BAM_LIST | tail -n 1`
BAM_OUTPUT=`head -n $SGE_TASK_ID $BAM_LIST | tail -n 1 | cut -d. -f1`

java -Xmx2g -jar ${MARK_DUPS} \
INPUT=${BAM_INPUT} \
OUTPUT=${NO_DUPS_DIR}/${BAM_OUTPUT}.sorted.nodup.bam \
METRICS_FILE=${NO_DUPS_DIR}/${BAM_OUTPUT}.marked_dup_metrics.txt \
REMOVE_DUPLICATES=true