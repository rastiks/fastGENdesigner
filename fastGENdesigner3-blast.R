# INSTALLATION
# install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')
# download.file("https://ftp.ncbi.nlm.nih.gov/blast/db/human_genome.00.tar.gz","human_genome.00.tar.gz",  mode='wb')
# download.file("https://ftp.ncbi.nlm.nih.gov/blast/db/human_genome.01.tar.gz","human_genome.01.tar.gz",  mode='wb')
# untar("human_genome.00.tar.gz", exdir="human_genome")
# untar("human_genome.01.tar.gz", exdir="human_genome")


find_problem <- function(primer, sub_cl, chromosomes,i){
  target <- sub_cl[sub_cl$QueryID == primer &  # name
                     grepl(paste("NC_0+", chromosomes[i,1], "[.]\\d+", sep=""), sub_cl$SubjectID) &  # chromosome  
                     (sub_cl$S.start >= chromosomes$V2[i]) & (sub_cl$S.end <= chromosomes$V3[i]) ,] # coordinates
  # find max score
  max_score <- target$Bits
  
  # filter max score -6
  final_query <- sub_cl[(sub_cl$QueryID == primer) & (sub_cl$Bits >= max_score -6),]
  control <- data.frame(matrix(target, nrow = nrow(final_query), ncol = ncol(target),byrow=TRUE))
  control <- control == final_query
  final_query <- final_query[!( rowSums(control) == ncol(final_query)),]
  #final_query <- final_query[!( rowSums(t(apply(final_query, 1, function(x) target == x))) == ncol(final_query)),]
  return(list(target,final_query))
}

# BLAST
blasting <- function(input_file, output_dir, input_type,blast_db) {
  
  # blast initialization
  bl <-blast(db=blast_db)
  seq_for_blast<-readDNAStringSet(paste(output_dir,"primers.fasta",sep="/"))
  cl <- predict(bl, seq_for_blast, BLAST_args = '-task blastn-short')
  
  # adding the output to the excel file
  wb <- loadWorkbook(paste(output_dir,"fastGENdesigner-output.xlsx", sep ="/"))
  addWorksheet(wb, sheetName ="BLAST"); writeData(wb, sheet = "BLAST", cl)
  saveWorkbook(wb,paste(output_dir,"fastGENdesigner-output.xlsx", sep ="/"), overwrite = TRUE)
  
  # input A, B  - rozdelit tabulku po genoch - zistit na akom je gen chromozome
  discarding <- c()

  #genes <- unique(sapply(strsplit(names(seq_for_blast), split= "_"), "[[",1))
  #for (my_gene in genes) { # ma zmysel toto robit
  chromosomes <- read.table(list.files(output_dir, pattern = "full_seqs.bed", full.names=T))
  chromosomes <- chromosomes[, -c(5,6)] 
  chromosomes$V1 <- gsub("chr", "", chromosomes$V1)
  # TO DO - uplne zakazat () v inpute?
  chromosomes$V4 <-  gsub("[(].[)]","", chromosomes$V4)
  cl$QueryID <- gsub("[(].[)]","", cl$QueryID)
  
  # FOR LOOP - pre kazdu sekvenciu
  problems_df <- data.frame()
  for (i in 1:nrow(chromosomes)) {
    sub_cl <- cl[grepl(chromosomes[i,4],cl$QueryID),]
    # for loop pre kazdy primer pair - p0 az p4 ...
    for (forward_primer in unique(grep("-F", sub_cl$QueryID, value=TRUE))) {
      result_forward <- find_problem(forward_primer, sub_cl, chromosomes,i)
      target <- result_forward[[1]]
      final_query <- result_forward[[2]]
      
      # je viac zhod - problem - skontrolovat aj reverse - ak je na rovnako zlom chromozome vyradime ho
      if (nrow(final_query) > 0) {
        result_reverse <- find_problem(gsub("-F","-R",unique(final_query$QueryID)), sub_cl, chromosomes,i)
        target_reverse <- result_reverse[[1]]
        final_query_reverse <- result_reverse[[2]]
        for (i_problem in 1:nrow(final_query)) {
          if (nrow(final_query_reverse[(final_query_reverse$SubjectID == final_query$SubjectID[i_problem]) & 
                                       !(grepl("NT_",final_query$SubjectID[i_problem])) & 
                                       ((abs(final_query_reverse$S.start - final_query$S.start[i_problem]))<10000) ,]) > 0)  { #  blizko seba
            reverse_problem <- final_query_reverse[(final_query_reverse$SubjectID == final_query$SubjectID[i_problem]) & # same offtarget chromosome
                                                     ((abs(final_query_reverse$S.start - final_query$S.start[i_problem]))<10000) ,]
            forward_problem <- final_query[i_problem,]
            problems_df <- rbind(problems_df,rbind(forward_problem, reverse_problem))
            discarding <- c(discarding,forward_primer, gsub("-F","-R",forward_primer))
            break
          }
        }
      }
    }
  }
  
  blasted <- seq_for_blast[!(names(seq_for_blast) %in% discarding)]
  writeXStringSet(blasted,paste(output_dir,"primers.fasta",sep="/"))
  
  cat("These primer pairs were excluded:\n")
  for (pair in unique(gsub("-[FR]", "", discarding))) cat(paste(pair, "\n"))
  
  # warning if all pairs were excluded
  before <- unique(gsub("_p\\d-[FR]", "", names(seq_for_blast)))
  after <- unique(gsub("_p\\d-[FR]","",names(blasted)))
  
  if (length(which(!(before %in% after))) > 0) cat(paste("!!! WARNING: THERE ARE NO PRIMER PAIRS LEFT FOR",before[!(before %in% after)], "!!!\n" ))
  missing_primers <- before[!(before %in% after)]
  return(problems_df)
  }
# adjusting excel
#wb <- loadWorkbook(paste(output_dir,"fastGENdesigner-output.xlsx", sep ="/"))
#my_data <- read.xlsx(wb,sheet = "Primer properties",colNames=FALSE)
#addStyle(wb, "Primer properties",rows=which(my_data$X1 %in% discarding), cols=1:5, style = createStyle(bgFill= "red"), gridExpand=TRUE, stack = FALSE)
#saveWorkbook(wb,paste(output_dir,"fastGENdesigner-output.xlsx", sep ="/"), overwrite = TRUE)

main <- function(input_file, output_dir, input_type, blast_db) {
  cat ("\nStarting Step3 - BLAST \nSearching for offtargets\n")
  problems_df <- blasting(input_file, output_dir, input_type, blast_db)
  return(problems_df)
}

# # ARGS
# args = commandArgs(trailingOnly=TRUE)
# 
# if (length(args)>0) {
#   for(i in 1:length(args)){
#     eval(parse(text=args[[i]]))
#   }
# } else {
#   input_file="/home/ppola/bva/fastgen_xpolak37/fastGENdesigner/inputs_outputs//inputs/input.txt"
#   output_dir="/home/ppola/bva/fastgen_xpolak37/fastGENdesigner/inputs_outputs"
#   input_type = "A"
#   blast_db = "/home/ppola/ncbi-blast-2.14.0+/blast/db/GCF_000001405.39_top_level"
# }

problems_df <- main(input_file, output_dir, input_type, blast_db)
