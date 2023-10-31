# fastGENdesigner
fastGENdesigner

Requirements:
- R v4.2
- Primer3: https://github.com/primer3-org/primer3
- Primer pooler: http://ssb22.user.srcf.net/pooler/pooler.tgz
- hg38.2bit: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/
- BlastDB: human_genome.00.tar.gz and human_genome.01.tar.gz from: https://ftp.ncbi.nlm.nih.gov/blast/db/ .. then extract them to one file

- !! PLEASE UPDATE THE PATHS IN CONFIG FILE !!

Updates:
- fastGENdesigner now includes:
	* Modul1 - Seq Selection
	* Modul2 - Primer3
	* Modul3 - BLAST 
	* Modul4 - Primer Pooler
	
Usage:
- R script available to run fastGENdesigner:
1. You can modify fastGENdesigner-input and run it in terminal:
```
Rscript fastGENdesigner.R
```
2. Or you can source fastGENdesigner.R in RStudio

### Possible inputs:  
1. GENE name - see example <em>input.txt</em>  
fastGENdesigner will design primers for **all exons** of the given gene.  

|gene|
|:----|
|H3-3A|

2. ROIs coordinates and gene name - see example <em>input_all.txt</em>  
fastGENdesigner will design primers for **ROIs**. 

|start|stop|gene|
|:----|:----|:----|
|38|94|PIK3CA|
|104|111|PIK3CA|
|118|118|PIK3CA|
|344|350|PIK3CA|
|391|391|PIK3CA|
|418|420|PIK3CA|

3. Chromosomal coordinates - see example <em>input_c.txt</em>.

|chrom|start|stop|name|
|:----|:----|:----|:----|
|chr1|119510718|119510754|BAT40(T)37|
|chr2|39309549|39309575|MONO-27(T)27|
|chr2|47414421|47414447|BAT26(A)27|
|chr2|95183614|95183636|NR24(T)23|


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

- Modul4:   
poolfiles: primer pooler output    
suggestion of the most suitable poolfile(s)

- fastGENdesigner-output.xlsx: Excel file with Primer properties    
and primer pairs distributions created by primer pooler



