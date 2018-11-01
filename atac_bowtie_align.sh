#!/bin/bash
    
#$ -o /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/mito_unalign.out
#$ -e /som/sborrego/201810_ATACSEQ_MB468_R8/qsub_reports/mito_unalign.err
#$ -q free64,som,asom 
#$ -pe openmp 32   
#$ -m beas            
#$ -ckpt blcr         

# set -euxo pipefail

module load bowtie2/2.2.7
module load samtools/1.0


GENOME_CHROMS=/som/sborrego/refs/hg38_chroms_only
CHR_MITO=/som/sborrego/refs/hg38_Mitochondial_DNA

EXP_DIR=/som/sborrego/201810_ATACSEQ_MB468_R8
TRIM_DATA=${EXP_DIR}/trim_data

bowtie2 -x ${CHR_MITO} --un unaln.fq -1 YOUR_READ1.fq -2 YOUR_READ2.fq | samtools view -Sb -> mito.bam

bowtie2 --threads 32 -x GENOME_CHROMS -1 unaln_1.fq -2 unaln_2.fq | samtools view -Sb -> YOUR_sample.bam