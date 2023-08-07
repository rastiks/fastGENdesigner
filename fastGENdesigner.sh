#!/bin/bash
# main script for fastGENdesigner tool

while [[ $# -gt 0 ]]; do
 case $1 in
    --input)
       input_file="$2"
       shift
       shift
       ;;
    --output)
       output_dir="$2"
       shift
       shift
       ;;
    --size_range)
       size_range="$2"
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

temp_file=$(mktemp)
exec > >(tee "$temp_file") 2>&1

echo "Starting fastGENdesigner with these arguments:"
# arguments control
if [ -z ${input_file+x} ]; then input_file=$(cat fastGENdesigner-input | grep -v '^#' | grep "input_file" | cut -d '=' -f 2); echo "Input file: ${input_file}"; fi
if [ -z ${output_dir+x} ]; then output_dir=$(cat fastGENdesigner-input | grep -v '^#' | grep "output_dir" | cut -d '=' -f 2); echo "Output dir: ${output_dir}"; fi
if [ -z ${size_range+x} ]; then size_range=$(cat fastGENdesigner-input | grep -v '^#' | grep "size_range" | cut -d '=' -f 2); echo "Size range: ${size_range}"; fi
if [ -z ${pools_num+x} ]
then 
	pools_num=$(cat fastGENdesigner-input | grep -v '^#' | grep "pools" | cut -d '=' -f 2)
	if [ "${#pools_num}" == "0" ]
	then
		echo "Number of pools: primer pooler suggestion"
	else
		echo "Number of pools: ${pools_num}"
	fi
fi


primer3=$(cat config | grep "primer3=" | cut -d '=' -f 2)

input_type=$(Rscript fastGENdesigner1-seq_selection.R "input_file='$input_file'" "output_folder='$output_dir'" "comment=''")
comment=$(Rscript fastGENdesigner2-primer3.R "input_file='$input_file'" "output_folder='$output_dir'" "size_range='$size_range'" "input_type='$input_type'" "primer3_path='$primer3'")

if [ "${#comment}" != "0" ]
then 
	echo "Resizing sequences with unsuccessful primer design"
	input_type=$(Rscript fastGENdesigner1-seq_selection.R "input_file='${output_dir}/input_resizing.txt'" "output_folder='$output_dir'" "comment='$comment'")
	comment=$(Rscript fastGENdesigner2-primer3.R "input_file='${output_dir}/input_resizing.txt'" "output_folder='$output_dir'" "size_range='$size_range'" "input_type='$input_type'" "primer3_path='$primer3'")
fi

cat "$temp_file" >> ${output_dir}/output.log
rm "$temp_file"
bash fastGENdesigner3-primerpooler.sh --input $input_file --output $output_dir --pools $pools_num


