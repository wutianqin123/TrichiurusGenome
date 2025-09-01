#!/bin/sh
#####################################################
starttime=$(date +"%Y-%m-%d %T")
echo "Start Time:" $starttime
#####################################################
cd `dirname $0`
#cd /mnt/disk2/Lab_Users/huyun/Thailand_Rhinogobius/assemble_results/nf
	
#mafft_aln.pl \
#--dna_unalignedd nf_selected \
#--dna_aligned nf_aligned \
#--cpu 6

mafft_aln.pl \
	--dna_unaligned merged_nf \
	--dna_aligned merged_nf_aligned \
	--cpu 20

#####################################################
endtime=$(date +"%Y-%m-%d %T")
echo "End Time:" $endtime
#####################################################
start=$(date --date="$starttime" +%s);
end=$(date --date="$endtime" +%s);
seconds=$((end-start))
hour=$(( $seconds/3600 ))
min=$(( ($seconds-${hour}*3600)/60 ))
sec=$(( $seconds-${hour}*3600-${min}*60 ))
echo "Total running time:" ${hour}:${min}:${sec}
#####################################################
exit 0
