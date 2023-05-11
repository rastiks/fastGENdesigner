#!/bin/bash

echo "Starting fastGENdesigner"
input_file="/home/ppola/bva/fastgen/fastGENdesigner/input.txt"
output_folder="/home/ppola/bva/fastgen/fastGENdesigner"

Rscript fastGENdesigner1-seq_selection.R "input_file='$input_file'" "output_folder='$output_folder'"
Rscript fastGENdesigner2-primer3.R "input_file='$input_file'" "output_folder='$output_folder'"