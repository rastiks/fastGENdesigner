#!/bin/bash

while [[ $# -gt 0 ]]; do
 case $1 in
    --input)
       input_file="$2"
       shift
       shift
       ;;
    --output)
       output_folder="$2"
       shift
       shift
       ;;
    --pools)
       pools_num="$2"
       shift
       shift
       ;;
    *)
     echo "Input not understood."
       exit 1
       ;;
 esac
done

if [ "${#input_file}" == "0" ]; then input_file=$(cat fastGENdesigner-input | grep -v '^#' | grep "input_file" | cut -d '=' -f 2); fi
if [ "${#output_folder}" == "0" ]; then output_folder=$(cat fastGENdesigner-input | grep -v '^#' | grep "output_dir" | cut -d '=' -f 2); fi
if [ "${#pools_num}" == "0" ]; then pools_num=$(cat fastGENdesigner-input | grep -v '^#' | grep "pools" | cut -d '=' -f 2); fi
primer_pooler=$(cat config | grep "primer_pooler" | cut -d '=' -f 2)

temp_file=$(mktemp)
exec > >(tee "$temp_file") 2>&1

echo "Starting primer pooler"
# NUMBER OF POOLS - SUGGESTION
if [ "${#pools_num}" == "0" ]
then
	$primer_pooler --suggest-pools --genome=fastGENdesigner_files/hg38.2bit $output_folder/primers.fasta 2> ${output_folder}/pooler_output.txt
	pools_num=$(tail -1 ${output_folder}/pooler_output.txt | tr -dc '0-9' )
	#pools_num=$((pools_num-1)) # zmensit o jedna - prediskutovat
	rm ${output_folder}/pooler_output.txt
	rm overlap-report*.txt
	echo "Number of suggested pools: $pools_num"
else
	echo "Number of given pools: $pools_num"
fi

rm -rf ${output_folder}/pooler_output
mkdir -p ${output_folder}/pooler_output

# PRIMER POOLER
#$primer_pooler --pools=$pools_num,1,$output_folder/pooler_output/poolfile --genome=fastGENdesigner_files/hg38.2bit --seedless $output_folder/primers.fasta 2> ${output_folder}/pooler_output/pooler_output.txt
$primer_pooler --pools=$pools_num,1,$output_folder/pooler_output/poolfile --genome=fastGENdesigner_files/hg38.2bit $output_folder/primers.fasta 2> ${output_folder}/pooler_output/pooler_output.txt
rm overlap-report*.txt 

echo "Poolfiles and pooler output created"
echo ""

# POOLS DISTRIBUTION
good_pools=$(Rscript pooler_summary.R "input_file='$input_file'" "output_folder='$output_folder'")

if [ "${#good_pools}" != "0" ]
then
	echo "Selected pools:  ${good_pools}. Initiating scoring process"
	echo "Primer pooler scoring" > ${output_folder}/pooler_output/poolfiles_score.txt
	
	temp_file1=$(mktemp)
	for pool in $good_pools; do
		echo $pool >> ${output_folder}/pooler_output/poolfiles_score.txt
		$primer_pooler ${output_folder}/pooler_output/${pool} --print-bonds=5 >> ${output_folder}/pooler_output/poolfiles_score.txt 2> ${temp_file1}; 
	done
	rm "$temp_file1"
	
	scores=$(Rscript pooler_scoring.R "input_file='$input_file'" "output_folder='$output_folder'")
	read -ra scores_array <<< "$scores"
	read -ra pools_array <<< "$good_pools"

	echo "Calculated scores (smaller is better):"

	for ((i = 0; i < ${#scores_array[@]};i++)); do
		echo "${pools_array[i]}: ${scores_array[i]}" ; 

	done

fi	
echo "Done"

cat "$temp_file" >> ${output_folder}/output.log
rm "$temp_file"
