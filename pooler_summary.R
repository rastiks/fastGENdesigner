
my_fun <- function(col){
  return(TRUE %in% (col>1))
}

pools_picking <- function(final) {
  #  DVOJKY NEMAME
  if (!(TRUE %in% (final > 1))) {
    sums <- apply(final,2, sum)
    ktore_skorovat <- names(sums)[sums == nrow(final)]
    cat(ktore_skorovat)
    
    # DVOJKY MAME - ALE STALE SA DA NIECO VYBRAT
  } else if (FALSE %in% apply(final,2,my_fun)) { 
    okay_pools <- !(apply(final,2,my_fun))
    ktore_skorovat <- colnames(final)[okay_pools]
    message("There are still overlaps in your pools, but you can use these poolfiles to put all primers together?:")
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
    for (pool in my_poolfiles) message(pool)
    cat(NULL)
  } else {
    message("I cannot pick the best pools. Check them on your own, or set higher pools number.")
    cat(NULL)
    
  }
}

pooler_summary <- function(output_folder){
  files=list.files(path = output_folder, pattern="poolfile\\d", full.names = T) 
  
  file <- 1
  for (poolfile in files) {
    my_fasta <- readDNAStringSet(poolfile)
    names_fasta <- names(my_fasta)
    names_fasta <- gsub("[(].[)]","", names_fasta)
    # TO DO - GENERALIZOVAT
    #primer_pairs <- unique(regmatches(names_fasta,regexpr("_\\d+_\\d+_p",names_fasta)))
    primer_pairs <- unique(regmatches(names_fasta,regexpr(".+.\\d_*\\d*_p",names_fasta)))
    #primer_pairs <- unique(regmatches(names_fasta,regexpr("TARGET\\d",names_fasta)))
    if (file == 1) final=data.frame("primer_pairs"= primer_pairs)
    file <- file + 1 
    counts <- c()
    for (pair in primer_pairs) {
      counts <- c(counts,length(grep(pair, names_fasta))/2)
    }
    
    df <- data.frame(
      primer_pairs = primer_pairs,
      x = counts,
      stringsAsFactors = FALSE
    )
    colnames(df)[ncol(df)] <- basename(poolfile)
    final <- merge(final,df,by = "primer_pairs", all=TRUE)
  }
  
  final[is.na(final)] <- 0
  
  write.xlsx(final,paste(output_folder,"poolfiles_dist.xlsx", sep ="/"))
  
  final <- final[,-1]
  return(final)
}

main <- function(output_folder){
  message("Creating primer pooler summary")
  message("Loading packages")
  suppressMessages(library(seqinr))
  suppressMessages(library(openxlsx)) 
  suppressMessages(library(Biostrings))
  final <- suppressWarnings(pooler_summary(output_folder))
  message("Done")
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