#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("EnsDb.Hsapiens.v86") 


saving_files <- function(d, my_rois_final, output_folder,my_start,my_end, my_rois_seq) {
  if (length(d)  == 1) name_file <- "exons"
  else name_file <- "ROIs"
  
  # BED file - only exons
  df <- data.frame(seqnames=seqnames(my_rois_final),
                   starts=start(my_rois_final)-1,
                   ends=end(my_rois_final),
                   names=names(my_rois_final),
                   scores=c(rep(".", length(my_rois_final))),
                   strands=strand(my_rois_final))
  
  write.table(df, file=paste(output_folder,paste(name_file,".bed",sep=""),sep="/"), quote=F, sep="\t", row.names=F, col.names=F, append=TRUE)
  
  # BED file - -50 exons +50
  df_complete_seq <- data.frame(seqnames=seqnames(my_rois_final),
                                starts=my_start, 
                                ends=my_end,
                                names=names(my_rois_final),
                                scores=c(rep(".", length(my_rois_final))),
                                strands=strand(my_rois_final))
  
  
  write.table(df_complete_seq, file=paste(output_folder,paste(name_file,"_full_seqs.bed",sep=""),sep="/"), quote=F, sep="\t", row.names=F, col.names=F, append=TRUE)
  
  # FASTA
  names(my_rois_seq) <- names(my_rois_final)
  writeXStringSet(my_rois_seq,file = paste(output_folder,paste(name_file,"_full_sequences.fasta",sep=""),sep="/"),format="fasta", append=TRUE)
  
}

width_adjusting <- function(my_rois) {
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
  return (my_rois_final)
}

inverse_rle <- function(rle_object) {
  result <- c()
  for (i in 1:length(rle_object@values)) {
    result <- c(result,rep(rle_object@values[i],rle_object@lengths[i]))
  }
  return (result)
}


seq_selection <- function(input_file,output_folder) {
  
  # ENSEMBL DATABASE
  edb <- EnsDb.Hsapiens.v86
  
  # INPUTS
  d <- read.delim(input_file, header=T) 

  # <------------------- my_code
  if (length(d) == 4) {
    message("Chromosomal locations were given.")
    
    # in case of 
    if (sum(colnames(d) == c("chrom","start","stop", "name")) != 4) {
      first_row <- colnames(d)
      colnames(d) <- c("chrom","start","stop", "name")
      d <- rbind(first_row, d)
    }
    # TO DO  - rozdelit kod do podfunkcii - kontrola vstupu, rozdelenie dlhsich usekov
    my_rois <- GRanges(
      seqnames = Rle(d$chrom,rep(1,nrow(d))),
      ranges = IRanges(d$start, end=d$stop, names = d$name),
      strand = Rle(strand("+"), nrow(d))
    )
    
    # VLOZIT DO FUNKCIE
    prolong <- 100
    my_rois_final <- width_adjusting(my_rois)
    my_chrom <- inverse_rle(my_rois_final@seqnames)
    my_chrom <- paste("chr",my_chrom,sep="")
    my_start <- my_rois_final@ranges@start - prolong # problem
    my_end <- my_rois_final@ranges@start + my_rois_final@ranges@width + prolong 
    my_rois_seq <- getSeq(Hsapiens,my_chrom,start=my_start,end=my_end)
    my_rois_seq <- DNAStringSet(my_rois_seq)
    
    saving_files(d, my_rois_final, output_folder,my_start,my_end, my_rois_seq)
  }
  
  else {
  # MANE select
  message("These genes were given:")
  for (my_gene in unique(d$gene)) message(my_gene)
  message("")
  mane_gff <- readGFF("MANE.GRCh38.v1.0.ensembl_genomic.gff.gz")
  for (my_gene in unique(d$gene)) {
    message(paste("Working on ", my_gene))
    message("Searching for MANE SELECT")# check if gene name is correct
    
    if (my_gene %in% mane_gff$gene_name) message("MANE SELECT found")
    else {
      #message("Cannot find MANE SELECT, please try different gene name.")
      stop("Cannot find MANE SELECT, please try different gene name")
    }
    
    # ID of MANE SELECT
    prot_id <- unique(mane_gff[mane_gff$gene_name == my_gene,"protein_id"])
    prot_id <- prot_id[!is.na(prot_id)]
    
    trans_id <- unique(mane_gff[mane_gff$gene_name == my_gene,"transcript_id"])
    trans_id <- trans_id[!is.na(trans_id)]
    
    if (length(prot_id) == 1) {
      prot_id <- unlist(strsplit(prot_id, split="[.]"))[1]
      trans_id <- unlist(strsplit(trans_id, split="[.]"))[1]
      message(paste("Trancript ID:",trans_id, sep=" "))
      message(paste("Protein ID:",prot_id, sep=" "))
      message("Please check the MANE SELECT here:")
      message(paste("https://www.ensembl.org/Homo_sapiens/Location/View?db=core;g=ENSG00000163041;r=1:226062716-226072019;t=",trans_id,sep=""))
      message(paste("https://www.lrg-sequence.org/search/?query=",my_gene,sep=""))
      message(paste("https://www.genenames.org/data/gene-symbol-report/#!/symbol/",my_gene,sep=""))
      message(paste("https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=",my_gene,sep=""))
      message("")
    }
    else{
      # TO DO??
      stop("Found multiple MANE SELECT??, please choose the right one manually from:")
    }
    
    # WE WANT ALL EXONS
    if (length(d)  == 1) {
      # get all exons of our target gene
      message(paste("Checking exons for", my_gene,sep = " "))
      my_rois <- exons(edb, filter= ~ protein_id == prot_id)
      
      # delete UTR
      message("Removing UTRs")
      #message("-------------------------------------------------------")
      message("")
      five_utr <- (fiveUTRsByTranscript(edb, filter= ~ protein_id == prot_id))@unlistData
      three_utr <- (threeUTRsByTranscript(edb, filter= ~ protein_id == prot_id))@unlistData
      
      # check if some exon is all UTR - necessary because exon cannot by of length 0 - we have to delete whole exon 
      # 5' UTR
      if (TRUE %in% ((my_rois[my_rois$exon_id %in% five_utr$exon_id])@ranges@width == five_utr@ranges@width)) {
        where_five <- which(((my_rois[my_rois$exon_id %in% five_utr$exon_id])@ranges@width == five_utr@ranges@width) == TRUE)
        where_roi <- which(my_rois$exon_id == five_utr[where_five]$exon_id)
        my_rois <- my_rois[-where_roi]
        five_utr <- five_utr[-where_five]
      }
      # 3' UTR
      if (TRUE %in% ((my_rois[my_rois$exon_id %in% three_utr$exon_id])@ranges@width == three_utr@ranges@width)) {
        where_three <- which(((my_rois[my_rois$exon_id %in% three_utr$exon_id])@ranges@width == three_utr@ranges@width) == TRUE)
        where_roi <- which(my_rois$exon_id == three_utr[where_three]$exon_id)
        my_rois <- my_rois[-where_roi]
        three_utr <- three_utr[-where_three]
      }
      
      new_ranges_five <- IRanges(start = (five_utr@ranges@start - 5) + (five_utr@ranges@width), width = ((my_rois[my_rois$exon_id == five_utr$exon_id])@ranges@width) - five_utr@ranges@width + 5)
      my_rois[names((my_rois[my_rois$exon_id == five_utr$exon_id])@ranges)]@ranges <- new_ranges_five
      
      new_ranges_three <- IRanges(start = ((my_rois[my_rois$exon_id == three_utr$exon_id])@ranges@start) , width = ((my_rois[my_rois$exon_id == three_utr$exon_id])@ranges@width) - three_utr@ranges@width + 5)
      my_rois[names((my_rois[my_rois$exon_id == three_utr$exon_id])@ranges)]@ranges <- new_ranges_three
      
      # NAME THE OUTPUT - genename_ex_number_refseqID
      my_exon_nums <- c()
      for (n in 1:length(my_rois)) my_exon_nums<- c(my_exon_nums,mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"exon_number"][grep(my_rois$exon_id[n],(mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"exon_id"]))])
      refseq <- regmatches(unique(unlist(mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"Dbxref"])),regexpr("NM_\\d+[.]\\d",unique(unlist(mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"Dbxref"]))))
      names(my_rois) <- paste(my_gene,"ex",my_exon_nums, refseq,sep="_")
    }
    
    else{ 
      sub_input <- d[d$gene == my_gene,]
      prt <- IRanges(start=sub_input$start, end=sub_input$stop, names=rep(prot_id,nrow(sub_input)))
      gnm <- proteinToGenome(prt, edb)
      my_rois = NULL 
      my_rois = gnm[[1]]
      if (length(gnm)>1) {
        for(i in 2:length(gnm)) {
          my_rois<-c(my_rois,gnm[[i]]) 
        }
      }
      
      #names(my_rois) <- paste(my_rois$exon_id,my_rois$protein_start,my_rois$protein_end,sep='_') # HOW TO NAME THIS???
      my_exon_nums <- c()
      for (n in 1:length(my_rois)) my_exon_nums<- c(my_exon_nums,mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"exon_number"][grep(my_rois$exon_id[n],(mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"exon_id"]))])
      refseq <- regmatches(unique(unlist(mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"Dbxref"])),regexpr("NM_\\d+[.]\\d",unique(unlist(mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"Dbxref"]))))
      names(my_rois) <- paste(my_gene,"ex",my_exon_nums, refseq, my_rois$protein_start,my_rois$protein_end, sep="_")
    }
    my_rois_final <- width_adjusting(my_rois)
    
    # get sequences for primer design -100...exon....+100 or -200...roi....+200
    if (length(d)  == 1) prolong <- 100
    else prolong <- 200
    
    my_chrom <- rep(levels(my_rois_final@seqnames@values),length(my_rois_final))
    my_chrom <- paste("chr",my_chrom,sep="")
    my_start <- my_rois_final@ranges@start - prolong # problem
    my_end <- my_rois_final@ranges@start + my_rois_final@ranges@width + prolong 
    my_rois_seq <- getSeq(Hsapiens,my_chrom,start=my_start,end=my_end)
    my_rois_seq <- DNAStringSet(my_rois_seq)
    
    saving_files(d, my_rois_final, output_folder,my_start,my_end, my_rois_seq)
  }
  }
}

main <- function(input_file, output_folder){
  
  message("Starting step1 with args")
  message(paste0("input_file:", input_file))
  message(paste0("output_folder:", output_folder))
  
  message("Loading packages")
  message("")
  suppressMessages(library(EnsDb.Hsapiens.v86))
  suppressMessages(library(ensembldb))
  suppressMessages(library(seqinr))
  suppressMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
  seq_selection(input_file, output_folder)
  message("Done")
  message("")
}

# ARGS
args = commandArgs(trailingOnly=TRUE)

for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

#input_file="/home/ppola/bva/fastgen_xpolak37/fastGENdesigner/inputs_outputs/input_c.txt"
#output_folder="/home/ppola/bva/fastgen_xpolak37/fastGENdesigner/inputs_outputs"

main(input_file, output_folder)


