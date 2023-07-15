# fastGENdesigner
fastGENdesigner

R 4.2 required !

Primer3 required ! please install !
https://github.com/primer3-org/primer3

BLAST 2.9.0 + required
(sudo apt install ncbi-blast+)

Updates:
- Modul1 and Modul2 can be run via terminal with:
```
input_file="path_to_input/input.txt"
output_folder="path_to_output_dir"
Rscript fastGENdesigner1-seq_selection.R "input_file='$input_file'" "output_folder='$output_folder'"
Rscript fastGENdesigner2-primer3.R "input_file='$input_file'" "output_folder='$output_folder'"
```

- Bash script available to run modul1 and modul2, change of input_file and output_file needed!  

```
bash fastGENdesigner.sh
```

- UTRs removed
- automatic selection of MANE SELECT
- URL generator

### Possible inputs:  
1. GENE name - see example <em>input.txt</em>  
fastGENdesigner will design primers for **all exons** of the given gene.  

2. ROIs coordinates and gene name - see example <em>input_all.txt</em>  
fastGENdesigner will design primers for **ROIs**. 

3. Chromosomal coordinates - see example <em>input_c.txt</em>.


### Possible outputs:
- Modul1:  
1. **EXONS**  
exons.bed: bed file for all exons  
exons_full_seqs.bed: bed file for prolonged exons  
exons_full_sequences.fasta: prolonged exons sequences  

2. **ROIs**  
ROIs.bed: bed file for ROIs  
ROIs_full_seqs.bed: bed file for prolonged ROIs  
ROIs_full_sequences.fasta: prolonged ROIs sequences  

- Modul2:  
primers.bed: bed file for designed primers  
primers.fasta: primers sequences

- Modul3:
poolfiles - primer pooler output




