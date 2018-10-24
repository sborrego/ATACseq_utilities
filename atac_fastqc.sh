#!/bin/bash

#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/fastqc/fastqc.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/fastqc/fastqc.err 
#$ -q free64,som,asom
#$ -pe openmp 8-16
#$ -m beas
#$ -ckpt blcr

set -euxo pipefail

module load blcr
module load fastqc/0.11.7

DATA_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/hts.igb.uci.edu/sborrego18102289

QC_OUT_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/fastqc
QC_HTML_DIR=${QC_OUT_DIR}/fastqc_html

mkdir -p ${QC_OUT_DIR}
mkdir -p ${QC_HTML_DIR}

for FILE in `find ${DATA_DIR} -name \*.gz`; do
    fastqc $FILE \
    --outdir ${QC_OUT_DIR}

    mv ${QC_OUT_DIR}/*.html ${QC_HTML_DIR}
done

tar -C ${QC_OUT_DIR} -czvf ${QC_HTML_DIR}.tar.gz ${QC_HTML_DIR} 