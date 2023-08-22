
#' call primer3 for a given set of DNAstringSet object
#'
#' @param seq DNAstring object, one DNA string for the given amplicon
#' @param size_range default: '51-170'
#' @param Tm melting temprature parameters default:c(54,59,63)
#' @param name name of the amplicon in chr_start_end format
#' @param primer3 primer3 location
#' @param therme.param thermodynamic parameters folder
#' @param settings text file for parameters
#' @author Altuna Akalin modified Arnaud Krebs' original function
#' @example
#' 

# !!! primer 3 configuration NEEDED !!! 
.callP3NreadOrg<-function(input_list = list(seq = "ACGT",ROI_start = 100, ROI_width=3,size_range='50-170',Tm=c(54,59,63),
                          TmDiff = 3, name="unnamed", 
                          PolyX=6, PrimerSize=c(18,20,23), PrimerGC=c(30,50,70), comment=""),
                          primer3="~/primer3/primer3-2.6.1/src/primer3_core",
                          thermo.param="~/primer3/primer3-2.6.1/src/primer3_config/",
                          settings="fastGENdesigner_files/primer3_settings.txt"){
  
  # make primer 3 temporary input file
  p3.input=tempfile()
  p3.output=tempfile()
  ROI_target=paste(input_list$ROI_start,",",input_list$ROI_width)
  
  write(
    paste( sprintf("SEQUENCE_ID=%s\n",input_list$name),
           sprintf("SEQUENCE_TEMPLATE=%s\n",as.character(input_list$seq)),
           "PRIMER_TASK=pick_detection_primers\n",
           "PRIMER_PICK_LEFT_PRIMER=1\n" ,
           "PRIMER_PICK_INTERNAL_OLIGO=0\n",
           "PRIMER_PICK_RIGHT_PRIMER=1\n"  ,
           "PRIMER_EXPLAIN_FLAG=1\n"  ,
           sprintf("PRIMER_PAIR_MAX_DIFF_TM=%s\n", input_list$TmDiff),
           sprintf("PRIMER_MIN_TM=%s\n" ,input_list$Tm[1]),
           sprintf("PRIMER_OPT_TM=%s\n" ,input_list$Tm[2]),
           sprintf("PRIMER_MAX_TM=%s\n" ,input_list$Tm[3]),
           sprintf("PRIMER_MAX_POLY_X=%s\n" ,input_list$PolyX),
           sprintf("PRIMER_MIN_SIZE=%s\n" ,input_list$PrimerSize[1]),
           sprintf("PRIMER_MAX_GC=%s\n" ,input_list$PrimerGC[3]),
           sprintf("PRIMER_PRODUCT_SIZE_RANGE=%s\n" ,input_list$size_range),
           sprintf("SEQUENCE_TARGET=%s\n",ROI_target),    # width ,ROI_width 
           sprintf("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s\n" ,thermo.param),
           "=",
           sep=''
    )
    ,
    p3.input
  )
  
  # CALLING PRIMER3
  try(system(
    paste(primer3 ,p3.input, "--p3_settings_file",settings, 
          ">", p3.output)
  ))
  
  #import and parse the output into a dataframe 
  out=read.delim(p3.output, sep='=', header=FALSE)
  unlink(c(p3.input,p3.output) ) # delete temp files
  returned.primers=as.numeric(as.vector(out[out[,1]=='PRIMER_PAIR_NUM_RETURNED',][,2]))
  
  # check result
  if (length(returned.primers)==0){
    # SAVING STATISTICS if no primer was found
    #write.csv(out, file=paste(output_folder, "statistics_primer3.txt", sep="/"),append=TRUE)
    
    if (nchar(input_list$comment)>0) message(paste(input_list$comment,"...","\u2717", sep = " "))  
    else message(paste(input_list$name, "\u2717", sep = " "))
    return(NA)
  }
  
  if ((returned.primers)==0){ 
    # SAVING STATISTICS if no primer was found
    #write.table(out, file=paste(output_folder, "statistics_primer3.txt", sep="/"),append=TRUE)
    
    if (nchar(input_list$comment)>0) message(paste(input_list$comment,"...","\u2717", sep = " "))
    else message(paste(input_list$name, "\u2717", sep = " "))
    return(NA)
  }
  
  # PRIMERS SUCCESFULLY DESIGNED
  if (returned.primers>0){
    if (nchar(input_list$comment)>0) message(paste(input_list$comment,"...","\u2713", sep = " "))  
    else message(paste(input_list$name, "\u2713", sep = " "))
    
    designed.primers=data.frame()
    for (i in seq(0,returned.primers-1,1)){
      #IMPORT SEQUENCES
      id=sprintf(  'PRIMER_LEFT_%i_SEQUENCE',i)
      PRIMER_LEFT_SEQUENCE=as.character(out[out[,1]==id,][,2])
      id=sprintf(  'PRIMER_RIGHT_%i_SEQUENCE',i)
      PRIMER_RIGHT_SEQUENCE=as.character(out[out[,1]==id,][,2])
      
      #IMPORT PRIMING POSITIONS
      id=sprintf(  'PRIMER_LEFT_%i',i)
      PRIMER_LEFT=as.numeric(unlist(strsplit(as.vector((out[out[,1]==id,][,2])),',')))
      #PRIMER_LEFT_LEN=as.numeric(unlist(strsplit(as.vector((out[out[,1]==id,][,2])),',')))
      id=sprintf(  'PRIMER_RIGHT_%i',i)
      PRIMER_RIGHT=as.numeric(unlist(strsplit(as.vector((out[out[,1]==id,][,2])),',')))
      #IMPORT Tm
      id=sprintf(  'PRIMER_LEFT_%i_TM',i)
      PRIMER_LEFT_TM=as.numeric(as.vector((out[out[,1]==id,][,2])),',')
      id=sprintf(  'PRIMER_RIGHT_%i_TM',i)
      PRIMER_RIGHT_TM=as.numeric(as.vector((out[out[,1]==id,][,2])),',')
      
      res=out[grep(i,out[,1]),]
      extra.inf=t(res)[2,,drop=FALSE]
      colnames(extra.inf)=sub( paste("_",i,sep=""),"",res[,1])
      extra.inf=extra.inf[,-c(4:9),drop=FALSE] # remove redundant columns
      extra.inf=apply(extra.inf,2,as.numeric)
      
      #Aggegate in a dataframe
      primer.info=data.frame(i,
                             PRIMER_LEFT_SEQUENCE,PRIMER_RIGHT_SEQUENCE,
                             PRIMER_LEFT_TM, PRIMER_RIGHT_TM,
                             PRIMER_LEFT_pos=PRIMER_LEFT[1],
                             PRIMER_LEFT_len=PRIMER_LEFT[2], 
                             PRIMER_RIGHT_pos=PRIMER_RIGHT[1], 
                             PRIMER_RIGHT_len=PRIMER_RIGHT[2],
                             t(data.frame(extra.inf)))
      
      rownames(primer.info)=NULL
      designed.primers=rbind(designed.primers, primer.info)
      
    }
  }
  return(designed.primers)
}

resizing_step <- function(resizing,primers_table, input_type, output_folder){
  # delete primers that was already created - example: _1 was OK but _2 was not, we have to delete _1
  resizing <- gsub("(_\\d)*$", "",resizing)
  for (i in 1:length(resizing)) {
    if (length(grep(resizing[i], primers_table$seq_name))>0) {
      primers_table <- primers_table[-(grep(resizing[i], primers_table$seq_name)),]
    }
  }
  
  new_genes <- unique(gsub("_ex_.*","",resizing))
  if (input_type =="A") {
    new_input <- data.frame("gene"=new_genes)
    exons <- regmatches(resizing, regexpr(".+ex(_\\d+)", resizing)) 
    }
    write.csv(new_input, paste(output_folder,"input_resizing.txt", sep="/"), row.names=FALSE, quote = FALSE)
    return(list(primers_table, exons))
    
  # TO DO : INPUT TYPE B AND C
    
}

merging_files <- function(output_folder){
  files <- list.files(path = output_folder, pattern="(ROIs)|(exons)", full.names = T) 
  
  # primers - if there is A and B option - choose the best - the one that has all parts designed, or the one with more primers designed
  primers <- readDNAStringSet(paste(output_folder, "primers.fasta", sep="/"))
  options <- names(primers)[grep("_[AB]\\d_", names(primers))]
  removing <- c()
  for (seq in unique(gsub("_[AB]\\d_p\\d-[FR]","", options))) {
    if ((TRUE %in% grepl(paste(seq,"A1", sep="_"),options)) & (TRUE %in% grepl(paste(seq,"A2", sep="_"),options))) removing <- c(removing,names(primers[grep(paste(seq,"B\\d_p\\d", sep ="_"), names(primers))]))    
  }
  primers <- primers[!names(primers) %in% removing]
  writeXStringSet(primers,paste(output_folder, "primers.fasta", sep="/"))
  
  removing <- unique(gsub("\\d_p\\d-[FR]","", removing))
  primers <- read.csv(paste(output_folder, "primers.bed", sep="/"), sep="\t", header=FALSE)
  primers <- primers[!(gsub("\\d_p\\d+$","",primers$V4) %in% removing),]
  write.table(primers,paste(output_folder, "primers.bed", sep="/"),quote=FALSE,row.names=FALSE, col.names = FALSE, sep="\t")
  
  # orig bed
  orig <- read.csv(files[grep("^.bed$",basename(gsub("exons|ROIs","",files)))], sep = "\t", header=FALSE)
  orig_resize <- read.csv(files[grep("^_resizing.bed$",basename(gsub("exons|ROIs","",files)))], sep = "\t", header=FALSE)
  resize_names <- unique(gsub("_\\d+$","",orig_resize$V4))
  bool_var <- !(gsub("_\\d+$","",orig$V4) %in% unique(gsub("_[AB]*\\d+$","",orig_resize$V4)))
  orig <- orig[bool_var,]
  orig <- rbind(orig,orig_resize)
  orig <- orig[!(substring(orig$V4,1,nchar(orig$V4)-1) %in% removing),] # remove all A or Bs
  write.table(orig,files[grep("^.bed$",basename(gsub("exons|ROIs","",files)))],quote=FALSE,row.names=FALSE, col.names = FALSE, sep="\t")
  unlink(files[grep("^_resizing.bed$",basename(gsub("exons|ROIs","",files)))])
  
  # full seqs bed
  orig <- read.csv(files[grep("^_full_seqs.bed$",basename(gsub("exons|ROIs","",files)))], sep = "\t", header=FALSE)
  orig_resize <- read.csv(files[grep("^_resizing_full_seqs.bed$",basename(gsub("exons|ROIs","",files)))], sep = "\t", header=FALSE)
  orig <- orig[bool_var,]
  orig <- rbind(orig,orig_resize)
  orig <- orig[!(substring(orig$V4,1,nchar(orig$V4)-1) %in% removing),] # remove all A or Bs
  write.table(orig,files[grep("^_full_seqs.bed$",basename(gsub("exons|ROIs","",files)))],quote=FALSE,row.names=FALSE, col.names = FALSE, sep="\t")
  unlink(files[grep("^_resizing_full_seqs.bed$",basename(gsub("exons|ROIs","",files)))])
  
  # fasta
  orig <- readDNAStringSet(files[grep("^_full_sequences.fasta$",basename(gsub("exons|ROIs","",files)))])
  orig_resize <-readDNAStringSet(files[grep("^_resizing_full_sequences.fasta$",basename(gsub("exons|ROIs","",files)))])
  orig <- orig[bool_var,]
  orig <- c(orig,orig_resize)
  orig <- orig[!(substring(names(orig),1,nchar(names(orig))-1) %in% removing)] # remove all A or Bs
  writeXStringSet(orig,files[grep("^_full_sequences.fasta$",basename(gsub("exons|ROIs","",files)))])
  unlink(files[grep("^_resizing_full_sequences.fasta$",basename(gsub("exons|ROIs","",files)))])
  
  # input
  unlink(paste(output_folder, "input_resizing.txt", sep ="/"))
  }

############################## primer_design - using external program PRIMER3    ###############

primer3caller <- function(input_file, output_folder, size_range, input_type, primer3_path){
  
  # exons have padding_number 100 (input A), other input types have padding_number
  # 200 and their name is ROIs
  
  if (input_type =="A") {
    padding_number <- 100
    name_file <- "exons"} 
  else {
    padding_number <- 200
    name_file <- "ROIs"}
  
  # SECOND TRY of primer design, files are named differently
  if (grepl("input_resizing.txt", input_file)) name_file <- paste(name_file, "resizing", sep="_")
  
  # full sequences - for primer design
  seq_complete <- readDNAStringSet(paste(output_folder,"/",name_file,"_full_sequences.fasta", sep=""))
  
  # bed file - without padding
  bed_file <- read.table(paste(output_folder,"/",name_file,".bed", sep = ""))
  res <- GRanges(seqnames = bed_file$V1,ranges = IRanges(start = bed_file$V2,
                                                         end = bed_file$V3,
                                                         names = bed_file$V4))
  
  # width of ROIs/Exons
  width <- res@ranges@width
  
  # variables initialization
  list_of_primer_seq <- list()
  df_primers_complete <-data.frame()
  primers = NULL
  
  # for loop through bed file
  primers_table <- data.frame()
  resizing <- c()
  for(i in 1:length(res)) {
    # primer 3 settings
    my_list <- list(seq = seq_complete[i], ROI_start = padding_number, ROI_width = width[i], 
                    size_range = size_range, Tm = c(57,59,62), TmDiff = 3, name = names(seq_complete)[i], 
                    PolyX = 4, PrimerSize=c(18,20,23), PrimerGC=c(30,50,70), comment="")
    
    # check if sequence is not already in the resizing vector
    if (length(grep(substring(my_list$name,1,nchar(my_list$name)-1), resizing))>0) next
    
    primers<-.callP3NreadOrg(input_list = my_list,
                             primer3=primer3_path,
                             thermo.param=paste(substring(primer3_path,1,regexpr("/src/",primer3_path)+4),"primer3_config",sep="")) 
    
    # changing parameters if primer design was unsuccessful
    
    n_try <- 1
    while (!(is.integer(nrow(primers))) && (n_try != 5)) {
      if (n_try ==1) {my_list$PolyX <- 6; my_list$comment <- "Increasing the PRIMER_MAX_POLY_X (4->6)"}
      else if (n_try==2) {my_list$TmDiff <- 5; my_list$comment <- "Increasing the PRIMER_PAIR_MAX_DIFF_TM (3->5)"}
      else if (n_try==3) {my_list$PrimerSize <- c(17,20,27) ; my_list$comment <- "Decreasing the PRIMER_MIN_SIZE (18->17) and Increasing the PRIMER_MAX_SIZE (23->27)"}
      else if (n_try==4) {my_list$Tm <- c(54,59,62) ; my_list$comment <- "Decreasing the PRIMER_MIN_TM (57->54)"}
      
      primers<-.callP3NreadOrg(input_list = my_list,
                               primer3=primer3_path,
                               thermo.param=paste(substring(primer3_path,1,regexpr("/src/",primer3_path)+4),"primer3_config",sep=""))
      n_try <- n_try +1
    }
    
    # primer design unsuccessful
    if (!(is.integer(nrow(primers)))) {
      message("Primer design unsuccessful.")
      resizing <- c(resizing,my_list$name)
    }
    
    message("")
    
    # primer design successful
    if (is.integer(nrow(primers))) {
      # saving primers
      reps <- nrow(primers) 
      scores=c(rep(".", reps))     
      strands_left=c(rep("+", reps))
      strands_right=c(rep("-", reps))
      
      # add left and right primers to df_primers_complete dataframe
      # why - padding_number??
      df_primer_LEFT <- data.frame(seqnames=res[i]@seqnames@values,  start=primers$PRIMER_LEFT_pos+res[i]@ranges@start-padding_number, end=primers$PRIMER_LEFT_pos-padding_number+res[i]@ranges@start+primers$PRIMER_LEFT_len, names = paste(names(res)[i], paste("p",primers$i,sep=""), sep="_"),scores=scores,strand=strands_left)
      df_primer_RIGHT <- data.frame(seqnames=res[i]@seqnames@values,  start=1+primers$PRIMER_RIGHT_pos+res[i]@ranges@start-padding_number-primers$PRIMER_RIGHT_len, end=primers$PRIMER_RIGHT_pos+res[i]@ranges@start-padding_number + 1 , names = paste(names(res)[i], paste("p",primers$i,sep=""), sep="_"),scores=scores,strand=strands_right)
      df_primer <- rbind(df_primer_LEFT,df_primer_RIGHT)
      df_primers_complete <- rbind(df_primers_complete,df_primer)     
      
      # sequences for primerpooler
      seq_name <- paste(df_primer_RIGHT$names,"-F", sep='')
      seq_name <- c(seq_name, paste(df_primer_LEFT$names,"-R", sep=''))
      
      list_of_primer_seq <- primers$PRIMER_LEFT_SEQUENCE
      list_of_primer_seq <- c(list_of_primer_seq, primers$PRIMER_RIGHT_SEQUENCE)
      
      # primer properties
      primers_table <- rbind(primers_table, 
                             data.frame(seq_name,list_of_primer_seq,
                                        nchar(list_of_primer_seq),
                                        c(df_primer_LEFT$start, df_primer_RIGHT$start),
                                        c(df_primer_LEFT$end, df_primer_RIGHT$end),
                                        c(primers$PRIMER_LEFT_TM,primers$PRIMER_RIGHT_TM),
                                        c(primers$PRIMER_LEFT_GC_PERCENT, primers$PRIMER_RIGHT_GC_PERCENT),
                                        c(primers$PRIMER_LEFT_HAIRPIN_TH, primers$PRIMER_RIGHT_HAIRPIN_TH),
                                        rep(primers$PRIMER_PAIR_PRODUCT_SIZE,2)))
      }
  } # end of the for loop - for all sequences
  
  # calling resizing part 
  if (length(resizing)>0 & !(grepl("input_resizing.txt", input_file))) { 
    # TO DO - druha podmienka prezatial musela zostat, do buducna volat resizing part niekolko krat
    # momentalne sa Modul1 spusti spatne len raz
    
    # changing output - removing designed primers for problematic sequences, returning problematic sequences in the for of comment
    resized_return <- resizing_step(resizing, primers_table, input_type, output_folder)
    primers_table <- resized_return[[1]]
    comment <- resized_return[[2]]
    cat(comment)
  }
  
  write.fasta(sequences = as.list(primers_table$list_of_primer_seq), names=primers_table$seq_name, file.out = paste(output_folder,'primers.fasta', sep="/"), open= "a")
  message("FASTA file created")
  write.table(df_primers_complete, file=paste(output_folder,"primers.bed", sep="/"), quote=F, sep="\t", row.names=F, col.names=F, append = TRUE)
  message("BED file created")
  
  # SAVING PRIMER PROPERTIES
  # if we are at resizing step - append data
  if (grepl("input_resizing.txt", input_file)) {
    wb <- loadWorkbook(paste(output_folder,"fastGENdesigner-output.xlsx", sep ="/"))
    my_data <- read.xlsx(wb,sheet = "Primer properties")
    writeData(wb, sheet = "Primer properties", primers_table, startRow=nrow(my_data)+2, colNames=FALSE)
  } else{
  
  # otherwise create new excel
  colnames(primers_table) <- c("Name", "Sequence","Primer Size", "Start", "End","Tm", "GC Content","Hairpin","Product Size")
  wb <- createWorkbook(); addWorksheet(wb, sheetName ="Primer properties"); writeData(wb, sheet = "Primer properties", primers_table)
  setColWidths(wb, "Primer properties", cols = 1:2, widths = "auto")
  }
  
  saveWorkbook(wb,paste(output_folder,"fastGENdesigner-output.xlsx", sep ="/"), overwrite = TRUE)
  message("Primers' characteristics saved to Excel file: fastGENdesigner-output.xlsx")
  
  # merging files
  if (grepl("input_resizing.txt", input_file)) {
    merging_files(output_folder)
    
  }
}

main <- function(input_file, output_folder, size_range, input_type, primer3_path, resizing=FALSE){
  message("Starting Step2 - primer3 -> Looking for suitable primers")
  suppressMessages(library(seqinr))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(Biostrings))
  suppressMessages(library(openxlsx)) 
  suppressWarnings(primer3caller(input_file, output_folder, size_range, input_type, primer3_path))
  message("")
}

# ARGS
args = commandArgs(trailingOnly=TRUE)

if (length(args)>0) {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
} else {
  # ONLY FOR TOOL DEVELOPMENT
  input_file="/home/ppola/bva/fastgen_xpolak37/fastGENdesigner/inputs_outputs/input.txt"
  output_folder="/home/ppola/bva/fastgen_xpolak37/fastGENdesigner/inputs_outputs"
  input_type="A"
  primer3_path="/home/ppola/primer3/primer3-2.6.1/src/primer3_core"
  size_range = "50-170"
}

main(input_file, output_folder, size_range, input_type, primer3_path)
