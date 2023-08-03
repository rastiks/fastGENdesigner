
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


# !!! primer 3 configuration NEEDED !!! <-----------------
.callP3NreadOrg<-function(seq,ROI_start = 100, ROI_width=3,size_range='50-170',Tm=c(54,59,63),name, 
                          primer3="/home/ppola/primer3/primer3-2.6.1/src/primer3_core",
                          thermo.param="/home/ppola/primer3/primer3-2.6.1/src/primer3_config/",
                          settings="/home/ppola/bva/fastgen_xpolak37/fastGENdesigner/primer3_settings.txt"){
  
  # make primer 3 input file
  p3.input=tempfile()
  p3.output=tempfile()
  ROI_target=paste(ROI_start,",",ROI_width)
  write(
    paste( sprintf("SEQUENCE_ID=%s\n",name  ),
           sprintf("SEQUENCE_TEMPLATE=%s\n",as.character(seq)),
           "PRIMER_TASK=pick_detection_primers\n",
           "PRIMER_PICK_LEFT_PRIMER=1\n" ,
           "PRIMER_PICK_INTERNAL_OLIGO=0\n",
           "PRIMER_PICK_RIGHT_PRIMER=1\n"  ,
           "PRIMER_EXPLAIN_FLAG=1\n"  ,
           "PRIMER_PAIR_MAX_DIFF_TM=3\n",
           sprintf("PRIMER_MIN_TM=%s\n" ,Tm[1]),
           sprintf("PRIMER_OPT_TM=%s\n" ,Tm[2]),
           sprintf("PRIMER_MAX_TM=%s\n" ,Tm[3]),
           sprintf("PRIMER_PRODUCT_SIZE_RANGE=%s\n" ,size_range),
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
  
  #import and parse the output into a dataframe named designed.primers
  out=read.delim(p3.output, sep='=', header=FALSE)
  
  unlink(c(p3.input,p3.output) ) # delete temp files
  returned.primers=as.numeric(as.vector(out[out[,1]=='PRIMER_PAIR_NUM_RETURNED',][,2]))
  if (length(returned.primers)==0){
    # SAVING STATISTICS if no primer was found
    write.csv(out, file=paste(output_folder, "statistics_primer3.txt", sep="/"),append=TRUE)
    return(NA)
    }
  if ((returned.primers)==0){ 
    # SAVING STATISTICS if no primer was found
    write.table(out, file=paste(output_folder, "statistics_primer3.txt", sep="/"),append=TRUE)
    return(NA)}
  
  if (returned.primers>0){
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
                             t(data.frame(extra.inf))
                             
      )
      rownames(primer.info)=NULL
      designed.primers=rbind(designed.primers, primer.info)
      
    }
    
    
  }
  return(designed.primers)
}


############################## primer_design - using external program PRIMER3    ###############

primer3caller <- function(input_file, path_to_files, size_range){
  d <- read.delim(input_file, header=T) 
  
  # if only gene in input -> exons wanted, else ROIs wanted
  if (length(d) == 1) {
    padding_number <- 100
    name_file <- "exons"} 
  else if (length(d) == 4) {
    padding_number <- 100
    name_file <- "ROIs"}
  else{
      padding_number <- 200
      name_file <- "ROIs"
    }
  
  # full padded sequences
  seq_complete <- readDNAStringSet(paste(path_to_files,"/",name_file,"_full_sequences.fasta", sep=""))
  
  # bed file - without padding
  bed_file <- read.table(paste(path_to_files,"/",name_file,".bed", sep = ""))
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
  for(i in 1:length(res)) {
    # calling PRIMER3 for padded sequence with parameters: ROI start (padding_number), ROI width, size range of wanted region, melting temperatures
    primers<-.callP3NreadOrg(seq=seq_complete[i],name = names(seq_complete)[i], ROI_start = padding_number, ROI_width=width[i],size_range = size_range,Tm=c(52,59,65))  
    
    # check if primers were found
    if (!(is.integer(nrow(primers)))) {
      message(paste("For sequence number",i, names(res[i]), "no primer was found. Trying size range: 50-300" ))
      # TO DO: pridat algoritmus na spatne volanie
      primers<-.callP3NreadOrg(seq=seq_complete[i],name = "test", ROI_start = padding_number, ROI_width=width[i],size_range = "50-300",Tm=c(52,59,65)) 
      if (is.integer(nrow(primers))) {
        message("Primers found")
      }
    }
    
    if (is.integer(nrow(primers))) {
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
      
      write.fasta(sequences = as.list(list_of_primer_seq), names = seq_name,file.out = paste(path_to_files,'primers.fasta', sep="/"),open = "a")
      
      # primer properties
      primers_table <- rbind(primers_table, 
                             data.frame(seq_name,list_of_primer_seq,
                                        c(primers$PRIMER_LEFT_TM,primers$PRIMER_RIGHT_TM),
                                        c(primers$PRIMER_LEFT_GC_PERCENT, primers$PRIMER_RIGHT_GC_PERCENT),
                                        rep(primers$PRIMER_PAIR_PRODUCT_SIZE,2)))
      }
  }
  write.table(df_primers_complete, file=paste(path_to_files,"primers.bed", sep="/"), quote=F, sep="\t", row.names=F, col.names=F)
  
  # SAVING PRIMER PROPERTIES
  colnames(primers_table) <- c("Name", "Sequence", "Tm", "GC Content", "Product Size")
  write.xlsx(primers_table,paste(path_to_files,"primers_properties.xlsx", sep ="/"))
}

main <- function(input_file, output_folder, size_range){
  message("Starting step2 - primer3")
  message("Loading packages")
  suppressMessages(library(seqinr))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(Biostrings))
  suppressMessages(library(openxlsx)) 
  suppressWarnings(primer3caller(input_file, output_folder, size_range))
  message("Done")
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
  input_file="/home/ppola/bva/fastgen_xpolak37/fastGENdesigner/inputs_outputs/input_c.txt"
  output_folder="/home/ppola/bva/fastgen_xpolak37/fastGENdesigner/inputs_outputs"
}

if (length(args) < 3) {
  size_range = "50-170"
}

main(input_file, output_folder, size_range)