#!/bin/bash
    
#$ -q free64,som,asom,pub64 
#$ -pe make 16
#$ -R y 
#$ -t 1-6 
#$ -m beas            
#$ -ckpt blcr         

set -euxo pipefail
SGE_TASK_ID=1

if [ $# -ne 1 ]; then
    echo "usage: data_dir"
    exit 1
fi

# Provide directory name for alignment ($1) on command line

module load bowtie2/2.2.7
module load samtools/1.0

# Location of Bowtie2 genom indexes
GENOME_CHROMS=/som/sborrego/refs/hg38_chroms_only/hg38_chr_only

EXP_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8
ALIGN_DIR=${EXP_DIR}/alignments/"$1"
UNALIGN_MITO_DIR=${ALIGN_DIR}/mitochondrial_unaligned
ALIGN_CHR_DIR=${ALIGN_DIR}/chr_only_alignments

mkdir -p ${ALIGN_CHR_DIR}

# Logs to keep track of things
RUNLOG=${ALIGN_DIR}/runlog_chrOnly_align_${SGE_TASK_ID}.txt
echo "Run by `whoami` on `date`" > ${RUNLOG}

FLAG=${ALIGN_DIR}/alignment_chrOnly_errors_${SGE_TASK_ID}.flagstat
echo "Run by `whoami` on `date`" > ${FLAG}

# File for lists of files to be processed
LIST_1=${UNALIGN_MITO_DIR}/read1_list.txt
LIST_2=${UNALIGN_MITO_DIR}/read2_list.txt

# Selecting file names using task ID
READ_1=${UNALIGN_MITO_DIR}/`head -n $SGE_TASK_ID $LIST_1 | tail -n 1`
READ_2=${UNALIGN_MITO_DIR}/`head -n $SGE_TASK_ID $LIST_2 | tail -n 1`

PREFIX=`basename ${READ_1} _READ1.fq.gz`
BAM_FILE=${ALIGN_CHR_DIR}/${PREFIX}.sorted.bam

# Message for file alignment
echo "Chromosome Only Alignment Summary for ${PREFIX}" >> ${RUNLOG}
echo "*** Aligning: ${PREFIX}" >> ${RUNLOG}

# Alignment of mitochondrial unaligned files. Alternative options: using -N 1, default is 0.
bowtie2 \
--threads 16 \
-x ${GENOME_CHROMS} \
-N 1 \
-1 ${READ_1} -2 ${READ_2} \
| samtools view -Sb - | samtools sort -@ 16 -o ${BAM_FILE} -O bam -T ${PREFIX} - \
2>> ${RUNLOG}

samtools index ${BAM_FILE}
samtools flagstat ${BAM_FILE} >> ${FLAG}
	
# bowtie2 -x ${CHR_MITO} --un-conc-gz ${UNALIGNED_MITO_FILE} -1 ${R1} -2 ${R2} | samtools view -Sb -> mito.bam
# bowtie2 --threads 32 -x GENOME_CHROMS -1 unaln_1.fq -2 unaln_2.fq | samtools view -Sb -> YOUR_sample.bam