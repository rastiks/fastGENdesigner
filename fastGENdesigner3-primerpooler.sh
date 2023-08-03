#!/bin/bash
input_file=$1
output_folder=$2

echo "Starting primer pooler"

# NUMBER OF POOLS - SUGGESTION
~/pooler/pooler --suggest-pools --genome=/home/ppola/pooler/hg38.2bit $output_folder/primers.fasta 2> ${output_folder}/pooler_output.txt
pools_num=$(tail -1 ${output_folder}/pooler_output.txt | tr -dc '0-9' )
pools_num=$((pools_num-1)) # zmensit o jedna - prediskutovat
rm ${output_folder}/pooler_output.txt
echo "Number of pools: $pools_num"

# PRIMER POOLER
~/pooler/pooler --pools=$pools_num,1,$output_folder/poolfile --genome=/home/ppola/pooler/hg38.2bit --seedless $output_folder/primers.fasta 2> ${output_folder}/pooler_output.txt
mv overlap-report-1.txt $output_folder

echo "You can find poolfiles and primerpooler output in: $output_folder"
echo ""

# POOLS DISTRIBUTION
good_pools=$(Rscript pooler_summary.R "input_file='$input_file'" "output_folder='$output_folder'")

if [ "${#good_pools}" != "0" ]
then
	echo "I have chosen ${good_pools}. Starting to score them."
	echo "Primer pooler scoring" > ${output_folder}/poolfiles_score.txt

	for pool in $good_pools; do
		echo $pool >> ${output_folder}/poolfiles_score.txt
		~/pooler/pooler ${output_folder}/${pool} --print-bonds=5 >> ${output_folder}/poolfiles_score.txt 2> ${output_folder}/pooler_output.txt; 
	done

	scores=$(Rscript pooler_scoring.R "input_file='$input_file'" "output_folder='$output_folder'")
	read -ra scores_array <<< "$scores"
	read -ra pools_array <<< "$good_pools"

	echo "Calculated scores (smaller is better):"

	for ((i = 0; i < ${#scores_array[@]};i++)); do
		echo "${pools_array[i]}: ${scores_array[i]}" ; 

	done

fi	
echo "Done"


