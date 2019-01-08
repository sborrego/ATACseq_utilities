#!/bin/bash
    
#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/encode_xCorrelate.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/encode_xCorrelate.err
#$ -q free64,som,asom,pub64 
#$ -pe make 16
#$ -R y 
#$ -t 1-6 
#$ -m as            
#$ -ckpt blcr

set -euxo pipefail

module load R/3.5.1
module load boost/1.63.0
module load samtools

# Loading phantompeakqualtools program
SPP=/som/sborrego/programs/phantompeakqualtools/run_spp.R

# Input files
TAG_ALIGN_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/analysis_encode/alignments/tagAlign
TAG_LIST=${TAG_ALIGN_DIR}/tag_list.txt

# Taskarray file name
TAG_INPUT=`head -n $SGE_TASK_ID $TAG_LIST | tail -n 1`
TAG_OUTPUT=`head -n $SGE_TASK_ID $TAG_LIST | tail -n 1 | cut -d. -f1 | cut -d/ -f 8 | cut -d_ -f 1`

# Output directory
PHANTOM_OUTPUT=/som/sborrego/201810_ATACSEQ_MB468_R8/analysis_encode/alignments/peak

# Using variable names atacseq pipeline
CC_SCORES_FILE="${PHANTOM_OUTPUT}/${TAG_OUTPUT}.cc.qc"
CC_PLOT_FILE="${PHANTOM_OUTPUT}/${TAG_OUTPUT}.cc.plot.pdf"

# CC_SCORE FILE format
# Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab> relPhantomPeakCoef <tab> QualityTag

Rscript ${SPP} \
-c=${TAG_INPUT} \
-p=8 \
-filtchr=chrM \
-savp=${CC_PLOT_FILE} \
-out=${CC_SCORES_FILE}

sed -r 's/,[^\t]+//g' ${CC_SCORES_FILE} > temp
mv temp ${CC_SCORES_FILE}


#Rscript run_spp.R -c=/som/sborrego/201810_ATACSEQ_MB468_R8/analysis_encode/alignments/tagAlign/4R109-L8-P2-CGTACTAG.filt.nodup.sample.15.MATE1.tagAlign -filtchr=chrM -savp=test.plot.pdf -out=test.cc.qc