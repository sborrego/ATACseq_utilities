#!/bin/bash

#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/fastqc_trim_clip.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/fastqc_trim_clip.err 
#$ -q free64,som,asom
#$ -pe openmp 6
#$ -m beas
#$ -ckpt blcr

set -euxo pipefail

module load blcr
module load fastqc/0.11.7

EXP_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8
TRIM_DATA=${EXP_DIR}/trim_data_clipped

QC_OUT_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/fastqc_trim_clipped
QC_HTML_DIR=${QC_OUT_DIR}/fastqc_trim_clipped_html

mkdir -p ${QC_OUT_DIR}
mkdir -p ${QC_HTML_DIR}

for FILE in `find ${TRIM_DATA} -name \*P.fq.gz`; do
    fastqc $FILE \
    --threads 6 \
    --outdir ${QC_OUT_DIR}

    mv ${QC_OUT_DIR}/*.html ${QC_HTML_DIR}
done

tar -C ${QC_OUT_DIR} -czvf ${QC_HTML_DIR}.tar.gz ${QC_HTML_DIR} 