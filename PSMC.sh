#!/bin/bash

###############################################
#
# PSMC.sh
#
# This shell is for automatic PSMC pipeline
#
# Please build bwa index and run "chmod 777 PSMC.sh" at first.
#
# Usage: ./PSMC.sh genome.fa num_of_thread R1.fq R2.fq run_folder sequence_depth psmc_PATH
#
# Written by Lu Liang
#
# 2021.2.4 at SHOU
#
###############################################

mkdir ${5}
#1.create run folder, build index and alignment
bwa mem -t ${2} -o ${5}/psmcOut.sam ${1} ${3} ${4} 
samtools view -b -@ ${2} -o ${5}/psmcOut.bam ${5}/psmcOut.sam
samtools sort -@ ${2} -o ${5}/psmcOut_sorted.bam ${5}/psmcOut.bam
echo "1.create run folder, build index and alignment processed"

#2.use samtools get diploid.fq.gz
samtools mpileup -C50 -uf ${1} ${5}/psmcOut_sorted.bam > ${5}/psmcOut.bcf
bcftools call -c ${5}/psmcOut.bcf > ${5}/psmcOut.vcf
depth2x=`expr ${6} \* 2`
vcfutils.pl vcf2fq -d ${7} -D ${depth2x} ${5}/psmcOut.vcf | gzip > ${5}/diploid.fq.gz
echo "2.use samtools get diploid.fq.gz processed"

#3.psmc
${7}/utils/fq2psmcfa -q20 ${5}/diploid.fq.gz > ${5}/diploid.psmcfa
${7}/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${5}/diploid.psmc ${5}/diploid.psmcfa
${7}/utils/psmc2history.pl ${5}/diploid.psmc | ${7}/utils/history2ms.pl > ${5}/ms-cmd.sh
${7}/utils/psmc_plot.pl ${5}/diploid ${5}/diploid.psmc
echo "3.psmc processed!"
