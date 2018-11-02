#!/bin/bash
    
#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/bowtie_mito_alignment_clip.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/bowtie_mito_alignment_clip.err
#$ -q free64,som,asom,pub64 
#$ -pe make 16
#$ -R y 
#$ -t 1-6 
#$ -m beas            
#$ -ckpt blcr         

set -euxo pipefail

if [ $# -ne 1 ]; then
    echo "usage: data_dir"
    exit 1
fi

# Provide directory name for alignment ($1) on command line

module load bowtie2/2.2.7
module load samtools/1.0

# Location of Bowtie2 genom indexes
GENOME_CHROMS=/som/sborrego/refs/hg38_chroms_only/hg38_chr_only
CHR_MITO=/som/sborrego/refs/hg38_Mitochondial_DNA/chrM

EXP_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8
DATA_DIR=${EXP_DIR}/"$1"

ALIGN_DIR=${EXP_DIR}/alignments/"$1"
ALIGN_MITO_DIR=${ALIGN_DIR}/mitochondrial_alignments
UNALIGN_MITO_DIR=${ALIGN_DIR}/mitochondrial_unaligned

mkdir -p ${ALIGN_DIR}
mkdir -p ${ALIGN_MITO_DIR}
mkdir -p ${UNALIGN_MITO_DIR}

RUNLOG=${ALIGN_DIR}/runlog_mito_align_${SGE_TASK_ID}.txt
echo "Run by `whoami` on `date`" >> ${RUNLOG}

ERRLOG=${ALIGN_DIR}/errlog_mito_align_${SGE_TASK_ID}.txt
echo "Run by `whoami` on `date`" >> ${ERRLOG}

LIST_1=${DATA_DIR}/read1_list.txt
LIST_2=${DATA_DIR}/read2_list.txt

READ_1=${DATA_DIR}/`head -n $SGE_TASK_ID $LIST_1 | tail -n 1`
READ_2=${DATA_DIR}/`head -n $SGE_TASK_ID $LIST_2 | tail -n 1`

# Generating file names and new paths
PREFIX=`basename ${READ_1} _1P.fq.gz`
UNALIGNED_MITO_FILE=${UNALIGN_MITO_DIR}/${PREFIX}_chrM_unaligned_READ%.fq.gz
ALIGNED_MITO_FILE=${ALIGN_MITO_DIR}/${PREFIX}_chrM_alignment.bam

# Message for file alignment
echo "*** Aligning: ${PREFIX}"
echo "Mitochondrial Alignment Summary for ${PREFIX}" >> ${RUNLOG}

# Mitochondrial DNA alignment with Bowtie2
bowtie2 \
--threads 16 \
-x ${CHR_MITO} \
--un-conc-gz ${UNALIGNED_MITO_FILE} \
-1 ${READ_1} -2 ${READ_2} \
| samtools view -Sb - > ${ALIGNED_MITO_FILE} \
2>> ${RUNLOG} \
1>> ${ERRLOG}


# bowtie2 -x ${CHR_MITO} --un-conc-gz ${UNALIGNED_MITO_FILE} -1 ${R1} -2 ${R2} | samtools view -Sb -> mito.bam
# bowtie2 --threads 32 -x GENOME_CHROMS -1 unaln_1.fq -2 unaln_2.fq | samtools view -Sb -> YOUR_sample.bam