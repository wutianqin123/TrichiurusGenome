#!/bin/sh
#####################################################
starttime=$(date +"%Y-%m-%d %T")
echo "Start Time:" $starttime
#####################################################
cd `dirname $0`


#cd /mnt/disk2/Lab_Users/huyun/Thailand_Rhinogobius/merged_nf_filtered_selected_aligned

concat_loci.pl \
--indir merged_nf_filtered_selected_aligned \
--outfile concat_goby

exit 0
