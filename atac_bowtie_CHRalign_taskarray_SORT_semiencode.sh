#!/bin/bash
    
#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/bowtie_chr_align_sort_filter_clip.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/bowtie_chr_align_sort_filter_clip.err
#$ -q free64,som,asom,pub64 
#$ -pe make 32
#$ -R y 
#$ -t 1-6
#$ -m beas            
#$ -ckpt blcr         

set -euxo pipefail

# Provide directory name for alignment ($1) on command line
module load bowtie2/2.2.7
module load samtools/1.0

# Location of Bowtie2 genom indexes
GENOME_CHROMS=/som/sborrego/refs/hg38_chroms_only/hg38_chr_only

EXP_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8

ALIGN_DIR=${EXP_DIR}/alignments/data_trim_clip_noSORT
UNALIGN_MITO_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/alignments/trim_data_clipped/mitochondrial_unaligned

ALIGN_DIR_2=${EXP_DIR}/analysis_encode/alignments
ALIGN_CHR_DIR=${ALIGN_DIR_2}/bowtie_output_noSort
ALIGN_SORT_DIR=${ALIGN_DIR_2}/sorted_alignments

RUNLOG_DIR=${ALIGN_DIR_2}/runlogs

mkdir -p ${ALIGN_CHR_DIR}
mkdir -p ${ALIGN_SORT_DIR}
mkdir -p ${RUNLOG_DIR}

# Logs to keep track of things
RUNLOG=${RUNLOG_DIR}/runlog_chrOnly_align_sort_filter_${SGE_TASK_ID}.txt
echo "Run by `whoami` on `date`" >> ${RUNLOG}

FLAG=${RUNLOG_DIR}/alignment_chrOnly_errors_sort_filter_${SGE_TASK_ID}.flagstat
echo "Run by `whoami` on `date`" >> ${FLAG}

# File for lists of files to be processed
LIST_1=${UNALIGN_MITO_DIR}/read1_list.txt
LIST_2=${UNALIGN_MITO_DIR}/read2_list.txt

# Selecting file names using task ID
# command to prepare read1 list: for FILE in *_READ1*; do echo `pwd`/${FILE} >> read1_list.txt; done
READ_1=`head -n $SGE_TASK_ID $LIST_1 | tail -n 1`
READ_2=`head -n $SGE_TASK_ID $LIST_2 | tail -n 1`

PREFIX=`basename ${READ_1} _READ1.fq.gz`
PREFIX2=`basename ${READ_1} _chrM_unaligned_READ1.fq.gz`

BAM_FILE=${ALIGN_CHR_DIR}/${PREFIX_2}.bam
SORTED_BAM_FILE=${ALIGN_SORT_DIR}/${PREFIX_2}.sorted.filtered.bam

# Message for file alignment
echo "Chromosome Only Sorted and Filtered Alignment Summary for ${PREFIX}" >> ${RUNLOG}
echo "*** Aligning: ${PREFIX}" >> ${RUNLOG}

# Alignment of mitochondrial unaligned files. 
# Alternative options: using -N 1, default is 0; --local, default is --end-to-end; -X2000, default is 500.
bowtie2 \
--threads 32 \
-x ${GENOME_CHROMS} \
--local \
-1 ${READ_1} -2 ${READ_2} \
-X2000 \
| samtools view -Sb - > ${BAM_FILE}
2>> ${RUNLOG}

samtools flagstat ${BAM_FILE} >> ${FLAG}

samtools view -F 1804 -f 2 -b ${BAM_FILE} \
| samtools sort -@ 32 -o ${SORTED_BAM_FILE} -O bam -T ${SORTED_BAM_FILE} - 
samtools index ${SORTED_BAM_FILE}
samtools flagstat ${SORTED_BAM_FILE} >> ${FLAG}
	
# bowtie2 -x ${CHR_MITO} --un-conc-gz ${UNALIGNED_MITO_FILE} -1 ${R1} -2 ${R2} | samtools view -Sb -> mito.bam
# bowtie2 --threads 32 -x GENOME_CHROMS -1 unaln_1.fq -2 unaln_2.fq | samtools view -Sb -> YOUR_sample.bam