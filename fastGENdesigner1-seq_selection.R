# INSTALLATION
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("EnsDb.Hsapiens.v86") 

input_check <- function(d) {
  if (length(d) == 1 & sum(toupper(colnames(d)) == toupper("gene"))==1) {
    input_type <- "A"
    colnames(d) <- "gene"
    
  }
  else if (length(d) == 4 & sum(toupper(colnames(d)) == toupper(c("chrom","start", "stop", "name")))==4) {
    input_type <- "C"
    colnames(d) <- c("chrom","start", "stop", "name")
    # BED -> R: bed indexes from 0, R from 1
    d$start <- d$start +1
  }
  else if (length(d) == 3 & sum(toupper(colnames(d)) == toupper(c("start","stop","gene"))) == 3) {
    input_type <- "B"
    colnames(d) <- c("start","stop","gene")
  }
  else stop("Unrecognized input") 
  
  my_return <- list("input_type"=input_type, "input"=d)
  return(my_return)
  
}

saving_files <- function(my_rois_final, output_dir, my_start, my_end, my_rois_seq, input_type, comment="") {
  # check the input type - because of different file names
  if (input_type=="A") name_file <- "exons"
  else name_file <- "ROIs"

  if (nchar(comment)>0) name_file <- paste(name_file, "resizing", sep="_")
  
  # BED file 
  df <- data.frame(seqnames=seqnames(my_rois_final),
                   starts=start(my_rois_final)-1, # because bed works with 0
                   ends=end(my_rois_final),
                   names=names(my_rois_final),
                   scores=c(rep(".", length(my_rois_final))),
                   strands=strand(my_rois_final))
  
  write.table(df, file=paste(output_dir,paste(name_file,".bed",sep=""),sep="/"), quote=F, sep="\t", row.names=F, col.names=F, append=TRUE)
  
  # BED file - full
  df_complete_seq <- data.frame(seqnames=seqnames(my_rois_final),
                                starts=my_start, 
                                ends=my_end,
                                names=names(my_rois_final),
                                scores=c(rep(".", length(my_rois_final))),
                                strands=strand(my_rois_final))
  
  
  write.table(df_complete_seq, file=paste(output_dir,paste(name_file,"_full_seqs.bed",sep=""),sep="/"), quote=F, sep="\t", row.names=F, col.names=F, append=TRUE)
  message("BED files created")
  
  # FASTA
  names(my_rois_seq) <- names(my_rois_final)
  writeXStringSet(my_rois_seq,file = paste(output_dir,paste(name_file,"_full_sequences.fasta",sep=""),sep="/"),format="fasta", append=TRUE)
  message("FASTA files created")
  }

width_adjusting <- function(my_rois, max_length=150, comment="") {
  # Function for splitting the longer exons than 150bp
  my_width <-my_rois@ranges@width
  longer_than <- which(my_width>max_length)
  problems <- my_rois[longer_than]
  
  if (length(problems) >0){
    my_rois <- my_rois[-longer_than]
    # message for user
    message(paste("These sequences are longer than ", max_length, "bp:", sep=""))
    message(paste(names(problems), collapse = "\n"))
    message("Starting width adjustment")
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
      split_number <- ceiling(problems@ranges@width[i]/max_length)
      
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

ratio_width_adjusting <- function(my_rois) {
    if (length(my_rois)>0) {
      my_starts <- c()
      my_names <- c()
      my_seqnames <- c()
      my_rois_id <- c()
      my_proteins_id <- c()
      my_strands <- c()
      my_starts <- c()
      my_widths <- c()
      # deal with longer exons/ROIs - split them so they are < 150bp
      for (i in 1:length(my_rois)){
        split_number <- ceiling(my_rois@ranges@width[i]/150)
        
        if (split_number == 2) { # lahsia situacia - rozdelime 60/40 a 40/60
        # exon names
        my_names <- c(my_names,paste(rep(names(my_rois[i]), factorial(split_number)*2),paste(rep(LETTERS[1:factorial(split_number)], each=2),as.character(seq(1,split_number)), sep=""),sep="_"))
        # chromosome
        my_seqnames <- c(my_seqnames,rep(levels(seqnames(my_rois[i])),factorial(split_number)*2))
        # exon id
        my_rois_id <- c(my_rois_id,rep(my_rois$exon_id[i], factorial(split_number)*2))
        # protein id
        my_proteins_id <- c(my_proteins_id,rep(my_rois$protein_id[i],factorial(split_number)*2))
        # strand
        my_strands <- c(my_strands,rep(as.character(my_rois[i]@strand@values),factorial(split_number)*2))
        
        # initial coordinates
        start <- my_rois@ranges@start[i]
        width <- my_rois@ranges@width[i]
        end <- start + width
        
        for (ratio in c(0.6,0.4)) {
          split_intervals <- c(start, start + round((end - start + 1) * ratio) - 1)
          
          my_starts <- c(my_starts,split_intervals)
          
          # new widths - diff 
          split_intervals <- c(split_intervals,end)
          my_widths <- c(my_widths, diff(split_intervals))
        }
  
        
      }
      
      # save all the data into Granges class
      gr <- GRanges(seqnames = my_seqnames, strand = my_strands,
                    ranges = IRanges(start = my_starts, width = my_widths),
                    protein_id = my_proteins_id, exon_id = my_rois_id)
      names(gr) <- my_names
      
      # merge the split problematic exons and nonproblematic ones
      }
    } else gr <- my_rois
  return(gr)
}

inverse_rle <- function(rle_object) {
  result <- c()
  for (i in 1:length(rle_object@values)) {
    result <- c(result,as.character(rep(rle_object@values[i],rle_object@lengths[i])))
  }
  return (result)
}

finalseq_creating <- function(input_type, output_dir, my_rois, comment){
  # width adjusting
  if (nchar(comment)>0) { 
    # vraciame sa z kroku 2 
    # predtym rozdelene neboli (kratsie ako 150bp - klasicky ich rozdelime)
    shorter_than <- (my_rois@ranges@width<=150)
    my_rois_final <- width_adjusting(my_rois[shorter_than],max_length = min(my_rois[shorter_than]@ranges@width) -1) # vyskusat funkcnost
    
    # predtym boli rozdelene - ostatne - rozdelime ich v inom pomere - TO DO 
    my_rois_final <- c(my_rois_final,ratio_width_adjusting(my_rois[!shorter_than]))
    
  } else my_rois_final <- width_adjusting(my_rois)
  
  # get sequences for primer design -100...exon....+100 or -200...roi....+200
  if (input_type  == "A") prolong <- 100
  else prolong <- 200
  
  my_chrom <- inverse_rle(my_rois_final@seqnames)
  my_chrom[!grepl("chr", my_chrom)] <- paste("chr", my_chrom[!grepl("chr", my_chrom)], sep="")
  my_start <- my_rois_final@ranges@start - prolong # problem
  my_end <- my_rois_final@ranges@start + my_rois_final@ranges@width + prolong 
  my_rois_seq <- getSeq(Hsapiens,my_chrom,start=my_start,end=my_end)
  my_rois_seq <- DNAStringSet(my_rois_seq)
  
  # saving created variables
  saving_files(my_rois_final, output_dir,my_start,my_end, my_rois_seq, input_type, comment)
  
}

seq_selection <- function(input_file,output_dir, comment) {
  # ----------------------------------------------------------------------------
  # MAIN function for MODULE1:
  # - input control
  # - MANE SELECT choice
  # - Sequence extraction
  # - UTR deletion
  # - Long sequences split
  # - Elongation for primer design
  
  # INPUTS:
  # input_file: path to input file
  # output_dir: path to output directory
  # comment: information about sequences that we are interested in, 
  # example: H3-3A_ex_3
  
  # OUTPUTS:
  # exons.bed / ROIS.bed: BED file of sequences that we are interested in
  # exons_full_seqs.bed / ROIs_full_seqs.bed: BED file of sequences that we are 
  # using for primer design
  # exons_full_sequences.fasta / ROIS_full_sequences.fasta: FASTA file of 
  # sequences for primer design
  
  # ----------------------------------------------------------------------------
  # ENSEMBL DATABASE
  edb <- EnsDb.Hsapiens.v86
  
  # INPUT
  d <- read.delim(input_file, header=T) 

  # input check - input d, output d and input_type
  input_check_return <- input_check(d)
  input_type <- input_check_return$input_type
  d <- input_check_return$input
  
  # INPUT C - chromosomal locations 
  if (input_type == "C") {
    message("Chromosomal locations were given")

    # Creating GRanges object
    my_rois <- GRanges(
      seqnames = Rle(d$chrom,rep(1,nrow(d))),
      ranges = IRanges(d$start, end=d$stop, names = d$name),
      strand = Rle(strand("+"), nrow(d))
    )
    
    # final step - elongation, saving
    finalseq_creating(input_type, output_dir, my_rois, comment)
  }
  
  else {
    # input A or B
      
    # MANE select
    message("These genes were given:")
    for (my_gene in unique(d$gene)) message(my_gene)
    message("")
    mane_gff <- readGFF("fastGENdesigner_files/MANE.GRCh38.v1.2.ensembl_genomic.gff.gz")
    
    # for loop in case of there are more genes than one given
    for (my_gene in unique(d$gene)) { 
      message(paste("Working on ", my_gene))
      message("Searching for MANE SELECT")# check if gene name is correct
      
      if (my_gene %in% mane_gff$gene_name) message("MANE SELECT found")
      else {
        stop("Cannot find MANE SELECT, please try different gene name")
      }
      
      # MANE SELECT ID
      prot_id <- unique(mane_gff[mane_gff$gene_name == my_gene,"protein_id"])
      prot_id <- prot_id[!is.na(prot_id)]
      
      trans_id <- unique(mane_gff[mane_gff$gene_name == my_gene,"transcript_id"])
      trans_id <- trans_id[!is.na(trans_id)]
      
      gene_id <-  unique(mane_gff[mane_gff$gene_name == my_gene,"gene_id"])
      gene_id <- gene_id[!is.na(gene_id)]
      
      if (length(prot_id) == 1) {
        # URL generator
        prot_id <- unlist(strsplit(prot_id, split="[.]"))[1]
        trans_id <- unlist(strsplit(trans_id, split="[.]"))[1]
        message(paste("Trancript ID:",trans_id, sep=" "))
        message(paste("Protein ID:",prot_id, sep=" "))
        message("Please check the MANE SELECT here:")
        message(paste("https://www.ensembl.org/Homo_sapiens/Location/View?db=core;g=",gene_id, sep=""))
        message(paste("https://www.lrg-sequence.org/search/?query=",my_gene,sep=""))
        message(paste("https://www.genenames.org/data/gene-symbol-report/#!/symbol/",my_gene,sep=""))
        message(paste("https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=",my_gene,sep=""))
        message("")
      }
      else{
        # TO DO??
        stop("Found multiple MANE SELECT??, please choose the right one manually from:")
      }
      
      # EXONS
      if (input_type == "A") {
        
        # get all exons 
        message(paste("Checking exons for", my_gene,sep = " "))
        my_rois <- exons(edb, filter= ~ protein_id == prot_id)
        
        # delete UTR
        message("Removing UTRs")
        message("")
        five_utr <- (fiveUTRsByTranscript(edb, filter= ~ protein_id == prot_id))@unlistData
        three_utr <- (threeUTRsByTranscript(edb, filter= ~ protein_id == prot_id))@unlistData
        
        # check if some exon is all UTR - necessary because exon cannot by of 
        # length 0 - we have to delete whole exon 
        # 5' UTR
        if (TRUE %in% ((my_rois[my_rois$exon_id %in% five_utr$exon_id])@ranges@width == five_utr@ranges@width)) {
          where_five <- which(((my_rois[my_rois$exon_id %in% five_utr$exon_id])@ranges@width == five_utr@ranges@width) == TRUE)
          where_roi <- which(my_rois$exon_id == five_utr[where_five]$exon_id)
          if (nchar(comment)==0) my_rois <- my_rois[-where_roi]
          five_utr <- five_utr[-where_five]
        }
        # 3' UTR
        if (TRUE %in% ((my_rois[my_rois$exon_id %in% three_utr$exon_id])@ranges@width == three_utr@ranges@width)) {
          where_three <- which(((my_rois[my_rois$exon_id %in% three_utr$exon_id])@ranges@width == three_utr@ranges@width) == TRUE)
          where_roi <- which(my_rois$exon_id == three_utr[where_three]$exon_id)
          my_rois <- my_rois[-where_roi]
          three_utr <- three_utr[-where_three]
        }
        
        # get only certain exons - resizing step from Module2
        if (nchar(comment)>0) { # attention - possible bug - I presume that exons are in order
          wanted  <- grep(my_gene,unlist(strsplit(comment, " ")), value=TRUE)
          wanted <- gsub(paste(my_gene,"_", sep=""),"", wanted)
          wanted <- as.integer(substring(wanted,4, nchar(wanted)))
          my_rois <- my_rois[wanted]
        }
        
        # removing UTR +-5bp
        new_ranges_five <- IRanges(start = (five_utr@ranges@start -5) + (five_utr@ranges@width), width = ((my_rois[my_rois$exon_id == five_utr$exon_id])@ranges@width) - five_utr@ranges@width + 5)
        my_rois[names((my_rois[my_rois$exon_id == five_utr$exon_id])@ranges)]@ranges <- new_ranges_five
        
        new_ranges_three <- IRanges(start = ((my_rois[my_rois$exon_id == three_utr$exon_id])@ranges@start) , width = ((my_rois[my_rois$exon_id == three_utr$exon_id])@ranges@width) - three_utr@ranges@width + 5)
        my_rois[names((my_rois[my_rois$exon_id == three_utr$exon_id])@ranges)]@ranges <- new_ranges_three
        
        # NAME OF THE OUTPUT - genename_ex_number_refseqID - possible bug - Exon ID do not match sometimes, why?
        my_exon_nums <- c()
        for (n in 1:length(my_rois)) {
          exon_num <- mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"exon_number"][grep(my_rois$exon_id[n],(mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"exon_id"]))]
          if (isEmpty(exon_num)) exon_num <- mane_gff[(mane_gff$gene_name == my_gene) & (mane_gff$type =="exon") & (mane_gff$start == my_rois@ranges@start[n]) & (mane_gff$end == (my_rois@ranges@start[n] + my_rois@ranges@width[n] -1)) ,"exon_number"]
          my_exon_nums <- c(my_exon_nums, exon_num)
          #my_exon_nums<- c(my_exon_nums,mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"exon_number"][grep(my_rois$exon_id[n],(mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"exon_id"]))])
        }  
        refseq <- regmatches(unique(unlist(mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"Dbxref"])),regexpr("NM_\\d+[.]\\d",unique(unlist(mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"Dbxref"]))))
        names(my_rois) <- paste(my_gene,"ex",my_exon_nums, refseq,sep="_")
        
      } else if (input_type =="B") {
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
        
        # NAME OF THE OUTPUT - genename_ex_number_refseqID - possible bug - Exon ID do not match sometimes, why?
        my_exon_nums <- c()
        for (n in 1:length(my_rois)) my_exon_nums<- c(my_exon_nums,mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"exon_number"][grep(my_rois$exon_id[n],(mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"exon_id"]))])
        refseq <- regmatches(unique(unlist(mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"Dbxref"])),regexpr("NM_\\d+[.]\\d",unique(unlist(mane_gff[(mane_gff$gene_name == my_gene) & mane_gff$type =="exon" ,"Dbxref"]))))
        names(my_rois) <- paste(my_gene,"ex",my_exon_nums, refseq, my_rois$protein_start,my_rois$protein_end, sep="_")
      }
      
      # COMMON STEP for input A and B
      finalseq_creating(input_type, output_dir, my_rois, comment)
      
    } # end of for loop
  }
  
  # returning input_type for other Modules
  return(input_type)
}

main <- function(input_file, output_dir, comment){
  message("")
  message("Starting Step1 - Seq Selection")
  if (nchar(comment)==0) input_type <- seq_selection(input_file, output_dir, comment)
  else input_type <- suppressMessages(seq_selection(input_file, output_dir, comment))
  message("")
  return(input_type)
}

# ARGS
# args = commandArgs(trailingOnly=TRUE)
# if (length(args) > 0) {
#   for(i in 1:length(args)){
#     eval(parse(text=args[[i]]))
#   }
# } else {
#   # ONLY FOR TOOL DEVELOPMENT
#   input_file="/home/ppola/bva/fastgen_xpolak37/fastGENdesigner/inputs_outputs/input_c.txt"
#   output_dir="/home/ppola/bva/fastgen_xpolak37/fastGENdesigner/inputs_outputs"
#   comment="H3-3A_ex_2 H3-3A_ex_4"
#   #input_file="/home/rastik/primer3/src/inputH3F3A.txt"
#   #output_dir="/home/rastik/primer3/src/outputs"
# }

input_type <- main(input_file, output_dir, comment)



