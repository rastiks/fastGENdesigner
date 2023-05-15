#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("EnsDb.Hsapiens.v86") 

seq_selection <- function(input_file,output_folder) {
  
  # ENSEMBL DATABASE
  edb <- EnsDb.Hsapiens.v86
  
  # INPUTS
  d <- read.delim(input_file, header=T) 

  # <------------------- my_code
  if (length(d)  == 1) {
    # get all exons of our target gene
    my_rois <- exons(edb, filter= ~ protein_id == d$gene[1])
    message(paste("Checking exons for", d$gene[1],sep = " "))
  }
  else{
    prt <- IRanges(start=d$start, end=d$stop, names=d$gene)
    gnm <- proteinToGenome(prt, edb)
    my_rois = NULL 
    my_rois = gnm[[1]]
    if (length(gnm)>1) {
      for(i in 2:length(gnm)) {
        my_rois<-c(my_rois,gnm[[i]]) 
      }
    }
    
    names(my_rois) <- paste(my_rois$exon_id,my_rois$protein_start,my_rois$protein_end,sep='_')
  }
  my_width <-my_rois@ranges@width
  longer_than_150 <- which(my_width>150)
  problems <- my_rois[longer_than_150]
  if (length(problems)>0) my_rois <- my_rois[-longer_than_150]
  
  if (length(problems) >0){
  # initialization
  my_starts <- c()
  my_names <- c()
  my_seqnames <- c()
  my_rois_id <- c()
  my_proteins_id <- c()
  my_strands <- c()
  my_starts <- c()
  my_widths <- c()
  
  # deal with longer exons/ROIs - split them so they are < 150bp
  for (i in 1:length(problems)){
    split_number <- ceiling(problems@ranges@width[i]/150)
    
    # exon names
    my_names <- c(my_names,paste(rep(names(problems[i]), split_number),as.character(seq(1,split_number)),sep="_"))
    # chromosome
    my_seqnames <- c(my_seqnames,rep(levels(seqnames(problems[i])),split_number))
    # exon id
    my_rois_id <- c(my_rois_id,rep(problems$exon_id[i], split_number))
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
                protein_id = my_proteins_id, exon_id = my_rois_id)
  names(gr) <- my_names
  
  # merge the split problematic exons and nonproblematic ones
  my_rois_final <- c(gr, my_rois)
  }
  else {
    my_rois_final <- my_rois
  }
  # get sequences for primer design -50...exon....+50 or -200...roi....+200
  if (length(d)  == 1) prolong <- 100
  else prolong <- 200
  
  my_chrom <- rep(levels(my_rois_final@seqnames@values),length(my_rois_final))
  my_chrom <- paste("chr",my_chrom,sep="")
  my_start <- my_rois_final@ranges@start - prolong # problem
  my_end <- my_rois_final@ranges@start + my_rois_final@ranges@width + prolong 
  my_rois_seq <- getSeq(Hsapiens,my_chrom,start=my_start,end=my_end)
  my_rois_seq <- DNAStringSet(my_rois_seq)
  
  if (length(d)  == 1) name_file <- "exons"
  else name_file <- "ROIs"
  
  
  # names=paste(gr_total$exon_id,gr_total$protein_start,gr_total$protein_end,sep='_')
  
  # BED file - only exons
  df <- data.frame(seqnames=seqnames(my_rois_final),
                   starts=start(my_rois_final)-1,
                   ends=end(my_rois_final),
                   names=names(my_rois_final),
                   scores=c(rep(".", length(my_rois_final))),
                   strands=strand(my_rois_final))
  
  write.table(df, file=paste(output_folder,paste(name_file,".bed",sep=""),sep="/"), quote=F, sep="\t", row.names=F, col.names=F)
  
  # BED file - -50 exons +50
  df_complete_seq <- data.frame(seqnames=seqnames(my_rois_final),
                                starts=my_start, 
                                ends=my_end,
                                names=names(my_rois_final),
                                scores=c(rep(".", length(my_rois_final))),
                                strands=strand(my_rois_final))
  
  
  write.table(df_complete_seq, file=paste(output_folder,paste(name_file,"_full_seqs.bed",sep=""),sep="/"), quote=F, sep="\t", row.names=F, col.names=F)
  
  # FASTA
  names(my_rois_seq) <- names(my_rois_final)
  writeXStringSet(my_rois_seq,file = paste(output_folder,paste(name_file,"_full_sequences.fasta",sep=""),sep="/"),format="fasta")
  # <------------------------
}

main <- function(input_file, output_folder){
  
  
  message("Starting step1 with args")
  message(paste0("input_file:", input_file))
  message(paste0("output_folder:", output_folder))
  
  message("Loading packages")
  suppressMessages(library(EnsDb.Hsapiens.v86))
  suppressMessages(library(ensembldb))
  suppressMessages(library(seqinr))
  suppressMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
  seq_selection(input_file, output_folder)
  message("Done")
}

# ARGS
args = commandArgs(trailingOnly=TRUE)

for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

#input_file="/home/ppola/bva/fastgen_xpolak37/fastGENdesigner/inputs_outputs/input_all.txt"
#output_folder="/home/ppola/bva/fastgen/fastGENdesigner/inputs_outputs"

main(input_file, output_folder)


