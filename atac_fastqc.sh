#!/bin/bash

#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/fastqc/fastqc.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/fastqc/fastqc.err 
#$ -q free64,som,asom
#$ -pe openmp 8-16
#$ -m beas
#$ -ckpt blcr

module load blcr
module load fastqc/0.11.7

# The directory where the data we want to analyze is located
DATA_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/hts.igb.uci.edu/sborrego18102289
# The directory where we want the result files to go
QC_OUT_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8/fastqc

# Making the result file directory
mkdir -p ${QC_OUT_DIR}

# Here we are performing a loop that will use each file in our data directory as input, "*" is a wild card symbol and in this context matches any file in the indicated directory
# Each file will be processed with the program "fastqc", "\" symbol indicates that more options for the program are on the next line 
# (--outdir) indicates the output directory for the result files
for FILE in `find ${DATA_DIR} -name \*.gz`; do
    fastqc $FILE \
    --outdir ${QC_OUT_DIR}
done