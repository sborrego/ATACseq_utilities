#!/bin/bash
    
#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/bowtie_alignment_clip.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/bowtie_alignment_clip.err
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

RUNLOG=${ALIGN_DIR}/runlog_align_${1}_${2}_${3}.txt
echo "Run by `whoami` on `date`" > ${RUNLOG}

FLAG=${ALIGN_DIR}/alignment_errors_${1}_${2}_${3}.flagstat
echo "Run by `whoami` on `date`" > ${FLAG}

LIST_1=${DATA_DIR}/read1_list.txt
LIST_2=${DATA_DIR}/read2_list.txt

READ_1=`head -n $SGE_TASK_ID $LIST_1 | tail -n 1`
READ_2=`head -n $SGE_TASK_ID $LIST_2 | tail -n 1`

for NUMBER in `seq $2 $3`; do
	for FILE in `find ${DATA_DIR} -name \*P"${NUMBER}"\*1P.fq.gz`; do
		PREFIX=`basename ${FILE} _1P.fq.gz`
		R1=${DATA_DIR}/"${PREFIX}"_1P.fq.gz
		R2=${DATA_DIR}/"${PREFIX}"_2P.fq.gz

		UNALIGNED_MITO_FILE=${UNALIGN_MITO_DIR}/${PREFIX}_chrM_unaligned_READ%.fq.gz
		ALIGNED_MITO_FILE=${ALIGN_MITO_DIR}/${PREFIX}_chrM_alignment.bam

		echo "*** Aligning: ${PREFIX}"
		echo "Mitochondrial Alignment Summary for ${PREFIX}" >> ${RUNLOG}

		bowtie2 \
		--threads 32 \
		-x ${CHR_MITO} \
		--un-conc-gz ${UNALIGNED_MITO_FILE} \
		-1 ${R1} -2 ${R2} \
		| samtools view -Sb - > ${ALIGNED_MITO_FILE} \
		2>> ${RUNLOG}
	done
done

ALIGN_CHR_DIR=${ALIGN_DIR}/chr_only_alignments

mkdir -p ${ALIGN_CHR_DIR}

for NUMBER in `seq $2 $3`; do
	for FILE in `find ${UNALIGN_MITO_DIR} -name \*P"${NUMBER}"\*1.fq.gz`; do
		
		PREFIX=`basename ${FILE} _READ1.fq.gz`
		PREFIX2=`basename ${FILE} _chrM_unaligned_READ1.fq.gz`
		R1=${UNALIGN_MITO_DIR}/"${PREFIX}"_READ1.fq.gz
		R2=${UNALIGN_MITO_DIR}/"${PREFIX}"_READ2.fq.gz

		BAM_FILE=${ALIGN_CHR_DIR}/${PREFIX2}.sorted.bam

		echo "*** Aligning: ${PREFIX}"
 		echo "Chromosome Only Alignment Summary for ${PREFIX}" >> ${RUNLOG}

		bowtie2 \
		--threads 32 \
		-x ${GENOME_CHROMS} \
		-N 1 \
		-1 ${R1} -2 ${R2} \
		| samtools view -Sb - | samtools sort -@ 32 -o ${BAM_FILE} -O bam -T ${PREFIX2} - 

		samtools index ${BAM_FILE}
		samtools flagstat ${BAM_FILE} >> ${FLAG}
	
	done
done

# bowtie2 -x ${CHR_MITO} --un-conc-gz ${UNALIGNED_MITO_FILE} -1 ${R1} -2 ${R2} | samtools view -Sb -> mito.bam
# bowtie2 --threads 32 -x GENOME_CHROMS -1 unaln_1.fq -2 unaln_2.fq | samtools view -Sb -> YOUR_sample.bam