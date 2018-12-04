#!/bin/bash
    
#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/bbmerge.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/bbmerge.err
#$ -q free64,som,asom,pub64 
#$ -pe make 16
#$ -R y 
#$ -t 1-12            
#$ -ckpt blcr  

# Including the option on SGE to email you will cause trouble if you have too many samples. 
# This may be more relevant to large samples like scRNAseq

set -euxo pipefail

if [ $# -ne 1 ]; then
    echo "usage: data_dir"
    exit 1
fi

# Provide directory name for alignment ($1) on command line

# Location of BBMerge
bbmerge=/data/apps/anaconda/3.7-5.3.0/bin/bbmerge.sh 

# Working directories
EXP_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8
DATA_DIR=${EXP_DIR}/"$1"
RESULTS=${EXP_DIR}/bbmerge_results

# Make non-existent files 
mkdir -p ${RESULTS}

# Logs to keep track of things
RUNLOG=${RESULTS}/runlog_chrOnly_align_${SGE_TASK_ID}.txt
echo "Run by `whoami` on `date`" >> ${RUNLOG}

FLAG=${RESULS}/alignment_chrOnly_errors_${SGE_TASK_ID}.flagstat
echo "Run by `whoami` on `date`" >> ${FLAG}

# File for lists of files to be processed
LIST_1=${DATA_DIR}/read1_list.txt
LIST_2=${DATA_DIR}/read2_list.txt

# Selecting file names using task ID
READ_1=${DATA_DIR}/`head -n $SGE_TASK_ID $LIST_1 | tail -n 1`
READ_2=${DATA_DIR}/`head -n $SGE_TASK_ID $LIST_2 | tail -n 1`

$bbmerge \
in1=$READ_1 in2=$READ_2 \
outa=${RESULTS}/adapters_${SGE_TASK_ID}.fa