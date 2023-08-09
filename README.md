# fastGENdesigner
fastGENdesigner

Requirements:
- R v4.2
- Primer3: https://github.com/primer3-org/primer3
- Primer pooler: http://ssb22.user.srcf.net/pooler/pooler.tgz
- hg38.2bit: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit

Updates:
- fastGENdesigner now includes:
	* Modul1 - Seq Selection
	* Modul2 - Primer3
	* Modul3 - BLAST (coming soon)
	* Modul4 - primerpooler.sh, pooler_summary.R, pooler_scoring.R
	
Usage:
- Bash script available to run fastGENdesigner
1. You can modify fastGENdesigner-input and run:
```
bash fastGENdesigner.sh
```
2. Or you can call fastGENdesigner like this:
```
bash fastGENdesigner.sh --input <input_file> --output <output_dir> --size_range 50-170 --pools 5
```

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
primer properties (Excel file)

- Modul4:
poolfiles - primer pooler output
primer pooler summary (Excel file)




