
my_fun <- function(col){
  return(TRUE %in% (col>1))
}

add_unique_worksheet <- function(wb, file, sheet_name) {
  existing_names <- getSheetNames(file)
  i <- 0
  new_name <- sheet_name
  while (new_name %in% existing_names) {
    i <- i+1
    new_name <- paste0(sheet_name, " (", i, ")")
  }
  addWorksheet(wb, sheetName = new_name)
  return(new_name)
}

combine_pools <- function(final) {
  tryCatch( 
    {
      okay_pools <- !(apply(final,2,my_fun))
      ktore_skorovat <- colnames(final)[okay_pools]
      #message("Some pools still have overlaps, but you can combine these pool files:")
      final <- final[,ktore_skorovat]
      my_poolfiles <- which.max(apply(final,2, sum))
      chybajuce <- which((final[,my_poolfiles] == 0)==TRUE)
      final <- final[,-my_poolfiles]
      
      while (length(chybajuce) > 0) {
        my_poolfile <- which.max(apply(final[chybajuce,],2,sum))
        my_poolfiles <- c(my_poolfiles,my_poolfile)
        chybajuce <- intersect(which((final[,names(my_poolfile)] == 0)==TRUE), chybajuce)
        final <- final[,-my_poolfile]
      }
      my_poolfiles <- names(my_poolfiles)
      message("Some pools still have overlaps, but you can combine these pool files:")
      for (pool in my_poolfiles) message(pool)
      cat(NULL)
    },
    error=function(e) {
      message("Unable to determine the best pools. Please review them manually or increase the number of pools")
    }
  )
}

pools_picking <- function(final) {
  if (!(TRUE %in% (final > 1))) {
    sums <- apply(final,2, sum)
    ktore_skorovat <- names(sums)[sums == nrow(final)]
    if (length(ktore_skorovat) == 0) {
      combine_pools(final)
    } else {
    cat(ktore_skorovat)
    }
    
  } else if (FALSE %in% apply(final,2,my_fun)) { 
      combine_pools(final)
    } else {
      message("Unable to determine the best pools. Please review them manually or increase the number of pools")
      cat(NULL)
    }
}

pooler_summary <- function(output_folder){
  files=list.files(path = paste(output_folder,"pooler_output", sep="/"), pattern="poolfile\\d", full.names = T) 
  
  file <- 1
  for (poolfile in files) {
    my_fasta <- readDNAStringSet(poolfile)
    names_fasta <- names(my_fasta)
    names_fasta <- gsub("[(].[)]","", names_fasta)
    # TO DO - GENERALIZOVAT
    #primer_pairs <- unique(regmatches(names_fasta,regexpr("_\\d+_\\d+_p",names_fasta)))
    primer_pairs <- unique(regmatches(names_fasta,regexpr(".+.\\d_*\\d*_p",names_fasta)))
    #primer_pairs <- unique(regmatches(names_fasta,regexpr("TARGET\\d",names_fasta)))
    if (file == 1) {
      final=data.frame("primer_pairs"= primer_pairs)
      final_names <- data.frame("primer_pairs"= primer_pairs)
    }
    file <- file + 1 
    counts <- c()
    pair_names <- c()
    for (pair in primer_pairs) {
      pair_names <- c(pair_names,
                      paste(unique(
                        substring(
                          regmatches(names_fasta[grep(pair, names_fasta)], 
                                     regexpr("_p\\d-[FR]",names_fasta[grep(pair, names_fasta)])),
                                      2,3)),collapse=" "))
      
      counts <- c(counts,length(grep(pair, names_fasta))/2)
      #counts <- c(counts,length(pair_names))
    }
    
    df <- data.frame(
      primer_pairs = primer_pairs,
      x = counts,
      stringsAsFactors = FALSE
    )
    
    df_names <- data.frame(
      primer_pairs = primer_pairs,
      x = pair_names,
      stringsAsFactors = FALSE
    )
    
    colnames(df)[ncol(df)] <- basename(poolfile)
    colnames(df_names)[ncol(df_names)] <- basename(poolfile)
    final <- merge(final,df,by = "primer_pairs", all=TRUE)
    final_names <- merge(final_names,df_names,by = "primer_pairs", all=TRUE)
  }
  
  final[is.na(final)] <- 0
  final_names[is.na(final_names)] <- "-"
  
  # adding dataframes into Excel file
  wb <- loadWorkbook(paste(output_folder,"fastGENdesigner-output.xlsx", sep ="/"))
  
  # adjusting primer properties - excluded primer pairs are marked as red
  before <- read.xlsx(wb, "Primer properties",colNames=TRUE)$Name
  after <- names(readDNAStringSet(paste(output_folder,"primers.fasta", sep ="/")))
  addStyle(wb, "Primer properties",rows=which(!(before %in% after))+1, cols=1:9, style = createStyle(bgFill= "red"), gridExpand=TRUE, stack = FALSE)
  sheet_name <- add_unique_worksheet(wb, file = paste(output_folder,"fastGENdesigner-output.xlsx", sep ="/"),
                                     sheet_name = paste("PrimerDistributionCounts-",length(files),sep ="")) 
  writeData(wb, sheet = sheet_name, final)
  setColWidths(wb, sheet_name, cols = 1, widths = "auto")
  
  sheet_name <- add_unique_worksheet(wb, file = paste(output_folder,"fastGENdesigner-output.xlsx", sep ="/"), 
                                     sheet_name = paste("PrimerDistributionNames-",length(files), sep =""))
  writeData(wb, sheet = sheet_name, final_names)
  setColWidths(wb,  sheet_name, cols = 1, widths = "auto")
  saveWorkbook(wb,paste(output_folder,"fastGENdesigner-output.xlsx", sep ="/"), overwrite = TRUE)
  
  
  message("Primers distribution saved to Excel file: fastGENdesigner-output.xlsx")
  final <- final[,-1]
  return(final)
}

main <- function(output_folder){
  message("Creating primer pooler summary")
  suppressMessages(library(seqinr))
  suppressMessages(library(openxlsx)) 
  suppressMessages(library(Biostrings))
  final <- suppressWarnings(pooler_summary(output_folder))
  message("Choosing the best pool(s).")
  pools_picking(final)
  message("")
}


args = commandArgs(trailingOnly=TRUE)

if (length(args)>0) {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
} else {
  # ONLY FOR TOOL DEVELOPMENT
  output_folder="/home/ppola/bva/fastgen_xpolak37/fastGENdesigner/inputs_outputs"
}

main(output_folder)