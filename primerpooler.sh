#!/bin/bash
input_file=$1
output_folder=$2

echo "Starting primer pooler"

pools_num=$(head -20 $output_folder/primers.fasta | grep "[-]F" | wc -l) # najst iny sposob

echo "Number of pools: $pools_num"

~/pooler/pooler --pools=5,1,$output_folder/poolfile --genome=/home/ppola/pooler/hg38.2bit $output_folder/primers.fasta

mv overlap-report-1.txt $output_folder
