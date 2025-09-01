#!/bin/sh
#####################################################
starttime=$(date +"%Y-%m-%d %T")
echo "Start Time:" $starttime
#####################################################
#indir可以敲全路径
cd `dirname $0`

#merge_loci.pl \
#--indir "/mnt/disk2/Lab_Users/huyun/Outgroup_for_Rhinogobius/assemble_results/nf /mnt/disk2/Lab_Users/huyun/Thailand_Rhinogobius/assemble_results/nf" \ 
#--outdir outdir \
#--min_seq 3

merge_loci.pl \
	--indir "/mnt/disk2/Lab_Users/huyun/Rhinogobius/Outgroup_for_Rhinogobius/assemble_results/nf /mnt/disk2/Lab_Users/huyun/Rhinogobius/Outgroup_for_Rhinogobius/ziling/assemble_results/nf /mnt/disk2/Lab_Users/huyun/Rhinogobius/Rhinogobius_China_MPE/2495/assemble_results/nf /mnt/disk2/Lab_Users/huyun/Rhinogobius/Rhinogobius_China_MPE/assemble_results/nf" \
	--outdir merged_nf \
	--min_seq 3

	
echo "Done!!!"

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
