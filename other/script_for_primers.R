suppressMessages(library(seqinr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Biostrings))

args = commandArgs(trailingOnly=TRUE)

for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

#input_file="/home/ppola/bva/fastgen_xpolak37/fastGENdesigner/inputs_outputs/orig_input.txt"
#output_folder="/home/ppola/bva/fastgen_xpolak37/fastGENdesigner/inputs_outputs"

my_primers <- readDNAStringSet(paste(output_folder,"primers.fasta", sep ="/"))

find_duplicates <- function(primers) {
  if (!(length(primers) == length(unique(primers)))) {
    duplicates <- names(primers[!names(primers) %in% names(unique(primers))])
    return(duplicates)
  }
}

for (i in seq(1,length(my_primers),by=10)) {# 5 because we have 5 primer pairs for each sequence/exon
  primers <- my_primers[i:(i+9)]
  
  # forward
  duplicates <- "x"
  pokus <- 0
  duplicates <- find_duplicates(primers[1:5])
  while (!is.null(duplicates)) {
    for (dup in duplicates) {
      #my_primers[dup] <- subseq(my_primers[dup],pokus+1,nchar(my_primers[dup])-pokus)
      #my_primers[dup] <- subseq(my_primers[dup],2,nchar(my_primers[dup])) # unique 4
      if (pokus > 0) my_primers[dup] <- subseq(my_primers[dup],1,nchar(my_primers[dup])-pokus)
      else my_primers[dup] <- subseq(my_primers[dup],2,nchar(my_primers[dup]))
    }
    
    primers <- my_primers[i:(i+9)]
    pokus <- pokus + 1
    duplicates <- find_duplicates(primers[1:5])
  }
  
  # reverse
  duplicates <- "x"
  pokus <- 0
  duplicates <- find_duplicates(primers[6:10])
  while (!is.null(duplicates)) {
    for (dup in duplicates) {
      #my_primers[dup] <- subseq(my_primers[dup],pokus+1,nchar(my_primers[dup])-pokus)
      #my_primers[dup] <- subseq(my_primers[dup],2,nchar(my_primers[dup])) # unique 4
      if (pokus > 0) my_primers[dup] <- subseq(my_primers[dup],1,nchar(my_primers[dup])-pokus)
      else my_primers[dup] <- subseq(my_primers[dup],2,nchar(my_primers[dup]))
    }
    primers <- my_primers[i:(i+9)]
    pokus <- pokus + 1
    duplicates <- find_duplicates(primers[6:10])
  }
}

write.fasta(as.list(my_primers),names(my_primers), paste(output_folder,"unique_primers.fasta", sep ="/"))