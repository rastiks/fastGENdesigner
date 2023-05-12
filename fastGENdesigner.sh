#!/bin/bash

echo "Starting fastGENdesigner"
input_file="/home/ppola/bva/fastgen/fastGENdesigner/inputs_outputs/input_all.txt"
output_folder="/home/ppola/bva/fastgen/fastGENdesigner/inputs_outputs"

Rscript fastGENdesigner1-seq_selection.R "input_file='$input_file'" "output_folder='$output_folder'"
Rscript fastGENdesigner2-primer3.R "input_file='$input_file'" "output_folder='$output_folder'"
