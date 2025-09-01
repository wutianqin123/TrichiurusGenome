#!/bin/sh
#####################################################
starttime=$(date +"%Y-%m-%d %T")
echo "Start Time:" $starttime
#####################################################
cd `dirname $0`

aligned=$1
ml=${aligned}_ml
partition=${aligned}_codon_partition
mkdir $ml

files=$(ls $aligned)

# parallel process
N=30
(
for file in $files
do
	((i=i%N));((i++==0)) && wait
	{
	raxmlHPC-SSE3 -f d -s $aligned/$file -m GTRGAMMA -n ${file%.fas}.tre -q $partition/${file%.fas}.txt -p 123456 -N 10 
	mv RAxML_bestTree.${file%.fas}.tre $ml/${file%.fas}.tre
	rm -rf RAxML_info.${file%.fas}*
	rm -rf RAxML_log.${file%.fas}*
	rm -rf RAxML_parsimonyTree.${file%.fas}*
	rm -rf RAxML_result.${file%.fas}*
	}&
done
)

echo "Done with RAxML!!!"

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
