#!/bin/sh
#####################################################
starttime=$(date +"%Y-%m-%d %T")
echo "Start Time:" $starttime
#####################################################
cd `dirname $0`

perl /mnt/genome8/Lab_Users/hjt_data/data_preparation/demultiplex_inline.pl \
--undemultiplexed  \
--demultiplexed /R_demultiplex \
--inline_index inlineindex.txt \
--index_pair indexpair.txt \
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
