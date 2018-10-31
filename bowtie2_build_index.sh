#!/bin/bash         

module load blcr
module load bowtie2/2.2.7

set -eoux pipefail

/som/sborrego/refs/hg38_chroms_only

REF_DIR=/som/sborrego/refs
CHROMFA=${REF_DIR}/hg38
CHR_DIR=${REF_DIR}/hg38_chroms_only

touch ${CHR_DIR}/hg38_chr_only.fa


## Buiding a Mitochondrial DNA Bowtie2 index
# tar -zxvf /som/sborrego/refs/hg38/hg38.chromFa.tar.gz ./chroms/chrM.fa
# bowtie2-build ./chroms/chrM.fa chrM

## Building a human genome, chromosome only Bowtie2 index
for NUM in `seq 1 22` "X" "Y" ; do
	FILE=./chroms/chr${NUM}.fa
	tar -xOzvf ${CHROMFA}/hg38.chromFa.tar.gz "${FILE}" >> ${CHR_DIR}/hg38_chr_only.fa
done

bowtie2-build ${CHR_DIR}/hg38_chr_only.fa hg38_chr_only