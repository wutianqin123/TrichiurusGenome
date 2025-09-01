#!/bin/sh
#####################################################
starttime=$(date +"%Y-%m-%d %T")
echo "Start Time:" $starttime
#####################################################
cd `dirname $0`

java -Xmx20g -jar ./astral.5.7.5.jar -i ./merged_nf_filtered_selected_aligned_ml/cat.tre -t 1 -o MPE_Rhinogobius_species_tree_constrained.tre

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
