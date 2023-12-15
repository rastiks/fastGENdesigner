# fastGENdesigner
list.of.packages <- c("seqinr", "Biostrings", "openxlsx","EnsDb.Hsapiens.v86",
                      "ensembldb","BSgenome.Hsapiens.UCSC.hg38","GenomicRanges",
                      "rBLAST")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if ( "rBLAST" %in% new.packages) {
  "Installing rBLAST - this can take a while ..\n"
  install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')
  new.packages <- new.packages[-(grep("rBLAST",new.packages))]
}

if(length(new.packages)>0) {
  cat("Installing packages - this can take a while ..\n")
  #install.packages(new.packages)
  BiocManager::install(new.packages)
}

time_start <- Sys.time()
if (interactive()) {
script_directory <- dirname(sys.frame(1)$ofile)
setwd(script_directory)
}
# ARGS
args = commandArgs(trailingOnly=TRUE)

if (length(args)>0) {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

inputs <- read.csv("fastGENdesigner-input", sep="=", comment.char = "#", header = FALSE,row.names =1 )
if (!(exists("input_file"))) input_file <- inputs["input_file",]
if (!(exists("output_dir"))) output_dir <- inputs["output_dir",]
if (!(exists("size_range"))) size_range <- inputs["size_range",]
if (!(exists("pools"))) num_pools <- inputs["pools",] else num_pools <- pools
if (!(exists("score_threshold"))) score_threshold <- as.integer(inputs["score_threshold",])
if (!(exists("dg_threshold"))) dg_threshold <- as.integer(inputs["dg_threshold",])
if (!(exists("attempts"))) num_attempts <- as.integer(inputs["attempts",]) else num_attempts <- as.integer(attempts)

zz <- file(paste(output_dir,"output.log", sep="/"), open = "wt")
sink(file = zz, type="output", split=TRUE)
#sink(file = zz, type="message", split=TRUE)
cat("Starting fastGENdesigner with these arguments:\n")
cat(paste("Input file:", input_file,"\n"))
cat(paste("Output dir:", output_dir),"\n")
cat(paste("Size range:", size_range),"\n")

if (is.na(num_pools)) {cat(paste("Number of pools:", "least as possible","\n"))} else {
  cat(paste("Number of pools:", num_pools,"\n"))}
if(is.na(score_threshold)) score_threshold <- 2
if(is.na(dg_threshold)) dg_threshold <- -2
if(is.na(num_attempts)) num_attempts <- 10

cat(paste("Attempts:", num_attempts),"\n")
cat(paste("Score threshold:", score_threshold),"\n")
cat(paste("dG threshold:", dg_threshold),"\n")

config <- read.csv("config", sep="=", comment.char = "#", header = FALSE,row.names =1)

cat("Loading packages - this can take a while ..\n")
suppressMessages(library(seqinr))
suppressMessages(library(Biostrings))
suppressMessages(library(openxlsx)) 
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(ensembldb))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
suppressMessages(library(GenomicRanges))
suppressMessages(library('rBLAST'))

# fastGENdesigner
comment <-""
source("fastGENdesigner1-seq_selection.R")

primer3_path <- config["primer3",]
source("fastGENdesigner2-primer3.R")
comment <- paste(comment, collapse = " ")

if (nchar(comment)>0) {
  input_file <- gsub(basename(input_file), "input_resizing.txt",input_file)
  cat("Resizing sequences with unsuccessful primer design\n")
  source("fastGENdesigner1-seq_selection.R")
  source("fastGENdesigner2-primer3.R")
}

blast_db <- config["blast_db",]
source("fastGENdesigner3-blast.R")

primer_pooler <- config["primer_pooler",]
source("fastGENdesigner4-primerpooler.R")

# output settings
# TO DO : time and date 
input_file <- inputs["input_file",]
d <- read.delim(input_file, header=T) 
time_end <- Sys.time()
duration <- paste(gsub("Time difference of ","",round(difftime(time_end,time_start),2)),"mins")
output_fastgen <- data.frame("Date"=time_start,
                             "Duration"=duration,
                             "Input File"=input_file,
                             "Output Directory"=output_dir,
                             "Number of attempts"=num_attempts,
                             "Number of pools"=num_pools,
                             "Score threshold"=score_threshold,
                             "dG threshold"=dg_threshold,
                             "Input Type"=input_type,
                             "Primer3"=primer3_path,
                             "Primer pooler"=primer_pooler,
                             "Blast DB"=blast_db,
                             "Final number of pools:"=number_of_parts)

if (is.null(combinations)){
  output_fastgen["Best pool(s)"] <- paste(success,collapse = " ")
} else {
  output_fastgen["Best combination"] <- paste(best_combination_list[[2]], collapse=" & ")
}

output_fastgen <- as.data.frame(t(output_fastgen))

if (is.null(best_combination_list)){ # 
  # create final fasta
  success <- unlist(success)
  final_file <- paste(output_dir,"pooler_output",success,sep="/")
  file.copy(final_file,output_dir)
  if (length(success)>1) renaming <- paste(output_dir,"/final_primers",regmatches(success,regexpr("\\d",success)),".fasta", sep = "")
  else renaming <- paste(output_dir,"final_primers.fasta", sep = "/")
  file.rename(paste(output_dir,success, sep="/"),renaming)
  final_primers_fasta <- list.files(output_dir,pattern = "final_primers")
  primers_names <- names(readDNAStringSet(paste(output_dir,final_primers_fasta,sep="/")))
} else{
  # create fasta from combinations <----  TO DO
  my_fasta <- readDNAStringSet(paste(output_dir,"primers.fasta",sep="/"))
  final_pick <- unlist(strsplit(best_combination_list[[2]],split=" [+] "))
  primers_names <- c()
  for (pick in final_pick){
    pick <- combinations[pick,]
    final_fasta <- my_fasta[gsub("-[FR]","",names(my_fasta)) %in% unlist(pick[,1:(ncol(pick)-2)])]
    primers_names <- c(primers_names,names(final_fasta))
    writeXStringSet(final_fasta,paste0(output_dir,"/","final_",rownames(pick),".fasta"))
  }
  combinations["Selected"] <- as.integer(rownames(combinations) %in% final_pick)
  best_combination_list[[1]]["Selected"] <- as.integer(rownames(best_combination_list[[1]]) %in% best_combination_list[[2]])
}
cat("FINAL FASTA FILE CREATED.\n")

# create bed
final_bed <- read.csv(paste(output_dir,"primers.bed",sep="/"),sep="\t",header=FALSE)
final_bed <- final_bed[final_bed[,4] %in% gsub("-[RF]","",primers_names),]
write.table(final_bed, file=paste(output_dir,"final_primers.bed",sep="/"), quote=F, sep="\t", row.names=F, col.names=F)
cat("FINAL BED FILE CREATED.\n")

excel_update(output_dir,list_dataframes[[1]],list_dataframes[[2]], success,combinations, best_combination_list,output_fastgen,d,final_bed[,4])

