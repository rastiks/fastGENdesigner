# fastGENdesigner
script_directory <- dirname(sys.frame(1)$ofile)
setwd(script_directory)

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
if (!(exists("pools"))) pools <- inputs["pools",]

#zz <- file(paste(output_dir,"output.log", sep="/"), open = "wt")
#sink(file = zz, type="output", split=TRUE)
#sink(file = zz, type="message", split=TRUE)
message("Starting fastGENdesigner with these arguments:")
message(paste("Input file:", input_file))
message(paste("Output dir:", output_dir))
message(paste("Size range:", size_range))

if (is.na(pools)) {message("Number of pools: number of final primer pairs")} else {
  message(paste("Number of pools:", pools))}

config <- read.csv("config", sep="=", comment.char = "#", header = FALSE,row.names =1)

message("Loading packages - this can take a while :)")
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
  message("Resizing sequences with unsuccessful primer design")
  source("fastGENdesigner1-seq_selection.R")
  source("fastGENdesigner2-primer3.R")
}

blast_db <- config["blast_db",]
source("fastGENdesigner3-blast.R")

primer_pooler <- config["primer_pooler",]
source("fastGENdesigner4-primerpooler.R")
