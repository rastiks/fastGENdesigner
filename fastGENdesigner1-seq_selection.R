#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("EnsDb.Hsapiens.v86") 

suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(ensembldb))
suppressMessages(library(seqinr))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg38"))

seq_selection <- function(input_file,output_folder) {
  
  # ENSEMBL DATABASE
  edb <- EnsDb.Hsapiens.v86
  
  # INPUTS
  d <- read.delim(input_file, header=T) 

  # <------------------- my_code
  if (length(d)  == 1) {
    # get all exons of our target gene
    my_exons <- exons(edb, filter= ~ protein_id == d$gene[1])
    # length of exons
    my_width <- my_exons@ranges@width
    # separate the problematic exons - longer than 150bp
    longer_than_150 <- which(my_width>150)
    problems <- my_exons[longer_than_150]
    my_exons <- my_exons[-longer_than_150]
  }
  
  # initialization
  my_starts <- c()
  my_names <- c()
  my_seqnames <- c()
  my_exons_id <- c()
  my_proteins_id <- c()
  my_strands <- c()
  my_starts <- c()
  my_widths <- c()
  
  # deal with longer exons - split them so they are < 150bp
  for (i in 1:length(problems)){
    split_number <- ceiling(problems@ranges@width[i]/150)
    
    # exon names
    my_names <- c(my_names,paste(rep(names(problems[i]), split_number),as.character(seq(1,split_number)),sep="_"))
    # chromosome
    my_seqnames <- c(my_seqnames,rep(levels(seqnames(problems[i])),split_number))
    # exon id
    my_exons_id <- c(my_exons_id,rep(problems$exon_id[i], split_number))
    # protein id
    my_proteins_id <- c(my_proteins_id,rep(problems$protein_id[i], split_number))
    # strand
    my_strands <- c(my_strands,rep(as.character(problems[i]@strand@values),split_number))
    
    # initial coordinates
    start <- problems@ranges@start[i]
    width <- problems@ranges@width[i]
    end <- start + width
    # spliting by width/split_number
    split_intervals <- round(seq(start,end,by=round(width/split_number)))
    
    # sometimes the last bp can be saved as new start - dont know why - can be deleted
    if (length(split_intervals) > split_number) split_intervals <- split_intervals[-length(split_intervals)]
    
    # new starts
    my_starts <- c(my_starts,split_intervals)
    
    # new widths - diff 
    split_intervals <- c(split_intervals,end)
    my_widths <- c(my_widths, diff(split_intervals))
  }
  
  # save all the data into Granges class
  gr <- GRanges(seqnames = my_seqnames, strand = my_strands,
                ranges = IRanges(start = my_starts, width = my_widths),
                protein_id = my_proteins_id, exon_id = my_exons_id)
  names(gr) <- my_names
  
  # merge the split problematic exons and nonproblematic ones
  my_exons_final <- c(gr, my_exons)
  
  # get sequences for primer design -50...exon....+50
  my_chrom <- rep(levels(my_exons_final@seqnames@values),length(my_exons_final))
  my_chrom <- paste("chr",my_chrom,sep="")
  my_start <- my_exons_final@ranges@start -50
  my_end <- my_start + my_exons_final@ranges@width + 50 
  my_exons_seq <- getSeq(Hsapiens,my_chrom,start=my_start,end=my_end)
  
  # BED file - only exons
  df <- data.frame(seqnames=seqnames(my_exons_final),
                   starts=start(my_exons_final)-1,
                   ends=end(my_exons_final),
                   names=names(my_exons_final),
                   scores=c(rep(".", length(my_exons_final))),
                   strands=strand(my_exons_final))
  
  write.table(df, file=paste(output_folder,"ROI_exons.bed",sep="/"), quote=F, sep="\t", row.names=F, col.names=F)
  
  # BED file - -50 exons +50
  df_complete_seq <- data.frame(seqnames=seqnames(my_exons_final),
                                starts=my_start, 
                                ends=my_end,
                                names=names(my_exons_final),
                                scores=c(rep(".", length(my_exons_final))),
                                strands=strand(my_exons_final))
  
  write.table(df_complete_seq, file=paste(output_folder,"ROI_exons_full_seqs.bed",sep="/"), quote=F, sep="\t", row.names=F, col.names=F)
  
  # FASTA
  names(my_exons_seq) <- names(my_exons_final)
  writeXStringSet(my_exons_seq,file = paste(output_folder,"exons_full_sequences.fasta",sep="/"),format="fasta")
  # <------------------------
}

main <- function(input_file, output_folder){
  message("Starting step1 with args")
  message(paste0("input_file:", input_file))
  message(paste0("output_folder:", output_folder))
  seq_selection(input_file, output_folder)
  message("Done")
}

# ARGS
args = commandArgs(trailingOnly=TRUE)

for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

main(input_file, output_folder)


