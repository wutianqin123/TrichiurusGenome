#!/bin/sh
#####################################################
starttime=$(date +"%Y-%m-%d %T")
echo "Start Time:" $starttime
#####################################################
cd `dirname $0`

#cd /mnt/disk2/Lab_Users/huyun/Thailand_Rhinogobius
#-s /mnt/disk2/Lab_Users/huyun/Thailand_Rhinogobius/merged_nf_filtered_selected_aligned \

/home/software/iqtree-1.7-beta9-Linux/bin/iqtree1.7 \
-s /mnt/disk2/Lab_Users/huyun/Rhinogobius/Rhinogobius_China_MPE/concat_goby.phy \
-bb \
1000 \
-nt AUTO \
-m MFP+MERGE \
-spp DNA_blocks1.txt \
-c 1 \
--prefix China_Rhinogobius

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
