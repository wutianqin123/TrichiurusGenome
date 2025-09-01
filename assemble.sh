#!/bin/sh
#####################################################
starttime=$(date +"%Y-%m-%d %T")
echo "Start Time:" $starttime
#####################################################
cd `dirname $0`
#cd /mnt/disk2/Lab_Users/huyun/Rhinogobius/Rhinogobius_China_MPE

assemble.pl \
--trimmed  \
--queryn  \
--queryp  \
--db "" \
--dbtype nucleo \
--ref_name  \
--outdir ./assemble_results \
--cpu 30

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
