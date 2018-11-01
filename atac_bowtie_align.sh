#!/bin/bash
    
#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/mito_unalign.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/mito_unalign.err
#$ -q free64,som,asom 
#$ -pe openmp 32   
#$ -m beas            
#$ -ckpt blcr         

set -euxo pipefail

# Provide directory name for alignment ($1) on command line

module load bowtie2/2.2.7
module load samtools/1.0

GENOME_CHROMS=/som/sborrego/refs/hg38_chroms_only
CHR_MITO=/som/sborrego/refs/hg38_Mitochondial_DNA

EXP_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8
DATA_DIR=${EXP_DIR}/"$1"

ALIGN_DIR=${EXP_DIR}/alignments/"$1"
ALIGN_MITO=${ALIGN_DIR}/mitochondrial_alignments
UNALIGN_MITO=${ALIGN_DIR}/mitchondrial_unaligned


mkdir -p ${ALIGN_DIR}


for NUMBER in `seq 1 6`; do
	for FILE in `find ${DATA_DIR} -name \*P"${NUMBER}"\*1P.fq.gz`; do
		PREFIX=`basename ${FILE} 1P.fq.gz`
		R1=${DATA_DIR}/"${PREFIX}"_1P.fq.gz
		R2=${DATA_DIR}/"${PREFIX}"_2P.fq.gz

		OUTPUT=${TRIM_DATA}/${PREFIX}_trimmed_clip.fq.gz


		echo ${R1}
		echo ${R2}
		echo ${OUTPUT}
		echo ${UNALIGN_MITO}/${PREFIX}
	done
done


bowtie2 -x ${CHR_MITO} --un-conc-gz ${UNALIGN_MITO}/${PREFIX} -1 YOUR_READ1.fq -2 YOUR_READ2.fq | samtools view -Sb -> mito.bam

bowtie2 --threads 32 -x GENOME_CHROMS -1 unaln_1.fq -2 unaln_2.fq | samtools view -Sb -> YOUR_sample.bam



mkdir -p ${QC_OUT_DIR}
mkdir -p ${QC_HTML_DIR}

for FILE in `find ${DATA_DIR} -name \*.gz`; do
    fastqc $FILE \
    --outdir ${QC_OUT_DIR}

    mv ${QC_OUT_DIR}/*.html ${QC_HTML_DIR}
done

tar -C ${QC_OUT_DIR} -czvf ${QC_HTML_DIR}.tar.gz ${QC_HTML_DIR} 