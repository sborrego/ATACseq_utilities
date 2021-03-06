#!/bin/bash
    
#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/trim_clip.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/trim_clip.err
#$ -q free64,som,asom 
#$ -pe openmp 32    
#$ -m beas            
#$ -ckpt blcr         

# set -euxo pipefail

module load blcr
TRIMMOMATIC_DIR=/data/apps/trimmomatic/0.35/trimmomatic-0.35.jar 
ADAPTERS=/data/apps/trimmomatic/0.35/adapters/NexteraPE-PE.fa

EXP_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8
DATA_DIR=${EXP_DIR}/hts.igb.uci.edu/sborrego18102289

TRIM_DATA=${EXP_DIR}/trim_data_clipped

mkdir -p ${TRIM_DATA}

RUNLOG=${TRIM_DATA}/runlog_trim_data.txt
echo "Run by `whoami` on `date`" > ${RUNLOG}

TRIMMER="ILLUMINACLIP:${ADAPTERS}:2:30:10 MINLEN:36"
# TRIMMER="CROP:93 HEADCROP:7 LEADING:3 TRAILING:1 SLIDINGWINDOW:4:15 MINLEN:36"

for NUMBER in `seq 1 6`; do
	for FILE in `find ${DATA_DIR} -name \*P"${NUMBER}"\*READ1-Sequences.txt.gz`; do
		PREFIX=`basename ${FILE} -READ1-Sequences.txt.gz`
		R1=${DATA_DIR}/"${PREFIX}"-READ1-Sequences.txt.gz
		R2=${DATA_DIR}/"${PREFIX}"-READ2-Sequences.txt.gz

		OUTPUT=${TRIM_DATA}/${PREFIX}_trimmed_clip.fq.gz

	    echo "*** ${PREFIX} Summary" >> $RUNLOG

	    java -jar ${TRIMMOMATIC_DIR} \
	    PE \
	    -threads 32 \
	    -phred33 \
	    ${R1} ${R2} \
	    -baseout ${OUTPUT} \
	    ${TRIMMER} \
	    2>> ${RUNLOG}

	done
done