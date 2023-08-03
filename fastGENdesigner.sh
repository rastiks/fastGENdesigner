#!/bin/bash

echo "Starting fastGENdesigner"
input_file="/home/ppola/bva/fastgen_xpolak37/fastGENdesigner/inputs_outputs/input_c.txt"
output_folder="/home/ppola/bva/fastgen_xpolak37/fastGENdesigner/inputs_outputs"
size_range="50-170"

Rscript fastGENdesigner1-seq_selection.R "input_file='$input_file'" "output_folder='$output_folder'"
Rscript fastGENdesigner2-primer3.R "input_file='$input_file'" "output_folder='$output_folder'" "size_range='$size_range'"
bash fastGENdesigner3-primerpooler.sh $input_file $output_folder
