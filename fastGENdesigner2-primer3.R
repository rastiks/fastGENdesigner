
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
#' primers=mapply(.callP3NreadOrg,seq=sample.seq,
#'                                   name=names(sample.seq),
#'                      MoreArgs=list(primer3=primer3,size_range=size_range,Tm=Tm,
#'                                    thermo.param=thermoParam,settings=settings)
#'                      ,SIMPLIFY = FALSE)
#' 
#' 

# pozor nutna konfiguracia dle cesty k primer3   !!!!!!!!!!!!!!!!!!!! 
.callP3NreadOrg<-function(seq,ROI_width=3,size_range='51-170',Tm=c(54,59,63),name,
                       primer3="/home/rastik/primer3/src/primer3_core",
                       thermo.param="/home/rastik/primer3/src/primer3_config/",
                       settings="/home/rastik/primer3/settings_files/default_settings.txt"){
  



  #print(excluded.regions)
  # make primer 3 input file
  p3.input=tempfile()
  p3.output=tempfile()
  ROI_target=paste("200,",ROI_width)
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
  #call primer 3 and store the output in a temporary file
  try(system(
    paste(primer3 ,p3.input, "-p3_settings_file",settings, 
                ">", p3.output)
  ))
  
  #import and parse the output into a dataframe named designed.primers
  out=read.delim(p3.output, sep='=', header=FALSE)
  
  unlink(c(p3.input,p3.output) ) # delete temp files
  returned.primers=as.numeric(as.vector(out[out[,1]=='PRIMER_PAIR_NUM_RETURNED',][,2]))
  if (length(returned.primers)==0){ warning('primers not detected for ',name,call. = FALSE);return(NA)}
  if ((returned.primers)==0){ warning('primers not detected for ',name,call. = FALSE);return(NA)}
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
      #print(primer.info)
    }
    
    
    
    #colnames(designed.primers)=c('PrimerID',
    #                             'Fwseq','Rvseq',
    #                             'FwTm','RvTm',
    #                             'FwPos','Fwlen',
    #                             'RvPos','Rvlen',
    #                             'fragLen' )
    
  }
  return(designed.primers)
}


#' convert primer list to genomic intervals
#' 
#' function returns a GRanges or data frame object from a list of primers designed
#' by \code{designPrimers} function after calculating genomic location of the amplicon
#' targeted by the primers.
#' 
#' @param primers a list of primers returned by \code{designPrimers} function
#' @param as.data.frame logical indicating if a data frame should be returned
#'        instead of \code{GRanges} object.
#'        
#' @examples
#'  data(bisPrimers)
#'  # remove data with no primers found
#'  bisPrimers=bisPrimers[!is.na(bisPrimers)]
#'  gr.pr=primers2ranges(bisPrimers) # convert primers to GRanges
#'          
#' @seealso \code{\link{filterPrimers}}, \code{\link{designPrimers}}
#'        
#' @export
#' @docType methods
primers2ranges<-function(primers,as.data.frame=FALSE){
  
  if(any(is.na(primers))){
    warning( "There are targets without primers\nfiltering those before conversion")
    primers=primers[ !is.na(primers) ]
  }
  df=do.call("rbind",primers) # get primers to a df
  locs=gsub("\\|\\.+","",rownames(df)) # get the coordinates from list ids
  temp=do.call("rbind",strsplit(locs,"_")) #
  
  start=as.numeric(temp[,2])
  chr=as.character(temp[,1])
  
  
  
  amp.start= start + as.numeric(df$PRIMER_LEFT_pos)  
  amp.end  = start + as.numeric(df$PRIMER_RIGHT_pos)
  res=data.frame(chr=chr,start=amp.start,end=amp.end,df)
  #saveRDS(res,file="/work2/gschub/altuna/projects/DMR_alignments/all.designed.primers.to.amps.rds")
  if(as.data.frame)
  {
    return(res)
  }
  gr=GRanges(seqnames=res[,1],ranges=IRanges(res[,2],res[,3]) )
  values(gr)=res[,-c(1,2,3)]
  gr
}


############################## primer_design - using external program PRIMER3    ###############
 library(seqinr)
 file.remove('primers.fasta')
 list_of_primers_df = list()
 list_of_primer_seq <- list()

 df_primers_complete <-data.frame()
 for(i in 1:length(gr_total)) 
     {
     primers<-.callP3NreadOrg(seq=seq_complete[i],name = "test", ROI_width=width[i],size_range = '51-170',Tm=c(52,59,65)) 
     new_member_of_list = list(primers)
     list_of_primers_df = c(list_of_primers_df,new_member_of_list)
     
      if (is.integer(nrow(primers))) { 
      reps <- nrow(primers) 
     
      scores=c(rep(".", reps))      # pozor nefunguje ked nenavrhne ziadny primer
      strands_left=c(rep("+", reps))
      strands_right=c(rep("-", reps))
     
     #obsahuje vzdy jeden ciel a par riadok bedu
     df_primer_LEFT <- data.frame(seqnames=gr_total[i]@seqnames@values,  start=primers$PRIMER_LEFT_pos+gr_total[i]@ranges@start-200, end=primers$PRIMER_LEFT_pos-200+gr_total[i]@ranges@start+primers$PRIMER_LEFT_len, names = paste(gr_total[i]$exon_id,gr_total[i]$protein_start,gr_total[i]$protein_end,primers$i,sep='_'),scores=scores,strand=strands_left)
     df_primer_RIGHT <- data.frame(seqnames=gr_total[i]@seqnames@values,  start=primers$PRIMER_RIGHT_pos+gr_total[i]@ranges@start-200-primers$PRIMER_RIGHT_len,end=primers$PRIMER_RIGHT_pos+gr_total[i]@ranges@start-200, names = paste(gr_total[i]$exon_id,gr_total[i]$protein_start,gr_total[i]$protein_end,primers$i,sep='_'),scores=scores,strand=strands_right)
     df_primer <- rbind(df_primer_LEFT,df_primer_RIGHT)
     df_primers_complete <- rbind(df_primers_complete,df_primer)     
      
     #export sekvencii do primerpoolera
     seq_name <- paste(gr_total[i]$exon_id,gr_total[i]$protein_start,gr_total[i]$protein_end,primers$i,"L", sep='_')
     seq_name <- c(seq_name, paste(gr_total[i]$exon_id,gr_total[i]$protein_start,gr_total[i]$protein_end,primers$i,"R", sep='_'))
     
     list_of_primer_seq <- primers$PRIMER_LEFT_SEQUENCE
     list_of_primer_seq <- c(list_of_primer_seq, primers$PRIMER_RIGHT_SEQUENCE)

     write.fasta(sequences = as.list(list_of_primer_seq), names = seq_name,file.out = 'primers.fasta',open = "a")
     }
   }

write.table(df_primers_complete, file="primers.bed", quote=F, sep="\t", row.names=F, col.names=F)
 
#####################################################################################################
