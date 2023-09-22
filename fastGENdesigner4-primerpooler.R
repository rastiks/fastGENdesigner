# STEP4 -PRIMER POOLER
message("Starting primer pooler")

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

excel_update <- function(output_dir,final=NULL,final_names=NULL, best_pools = NULL, combinations=NULL) {
  files <- list.files(path = paste(output_dir,"pooler_output", sep="/"), pattern="poolfile\\d", full.names = T) 
  # adding dataframes into Excel file
  wb <- loadWorkbook(paste(output_dir,"fastGENdesigner-output.xlsx", sep ="/"))
  
  # adjusting primer properties - excluded primer pairs are marked as red
  before <- read.xlsx(wb, "Primer properties",colNames=TRUE)$Name
  after <- names(readDNAStringSet(paste(output_dir,"primers.fasta", sep ="/")))
  addStyle(wb, "Primer properties",rows=which(!(before %in% after))+1, cols=1:9, style = createStyle(bgFill= "red"), gridExpand=TRUE, stack = FALSE)
  if (!(is.null(final))) {
    sheet_name <- add_unique_worksheet(wb, file = paste(output_dir,"fastGENdesigner-output.xlsx", sep ="/"),
                                       sheet_name = paste("PrimerDistributionCounts-",length(files),sep ="")) 
    writeData(wb, sheet = sheet_name, final,rowNames=TRUE)
    setColWidths(wb, sheet_name, cols = 1, widths = "auto")
    
    addStyle(wb, sheet_name,rows=1:(nrow(final)+1), cols=which(colnames(final) %in% best_pools) + 1,
             style = createStyle(bgFill= "green"), gridExpand=TRUE, stack = FALSE)
    
    sheet_name <- add_unique_worksheet(wb, file = paste(output_dir,"fastGENdesigner-output.xlsx", sep ="/"), 
                                       sheet_name = paste("PrimerDistributionNames-",length(files), sep =""))
    writeData(wb, sheet = sheet_name, final_names)
    setColWidths(wb,  sheet_name, cols = 1, widths = "auto")
  }
  
  if (!(is.null(combinations)) & is.data.frame(combinations)) {
    sheet_name <- add_unique_worksheet(wb, file = paste(output_dir,"fastGENdesigner-output.xlsx", sep ="/"),
                                       sheet_name = "Combinations to score") 
    writeData(wb, sheet = sheet_name, combinations,rowNames=TRUE)
    setColWidths(wb,  sheet_name, cols = 1:(ncol(combinations)+1), widths = "auto")
  } else{
    message("Primers distribution saved to Excel file: fastGENdesigner-output.xlsx")
  }
  
  
  saveWorkbook(wb,paste(output_dir,"fastGENdesigner-output.xlsx", sep ="/"), overwrite = TRUE)
  
}

pooler_summary <- function(output_dir){
  files=list.files(path = paste(output_dir,"pooler_output", sep="/"), pattern="poolfile\\d", full.names = T) 
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
  
  rownames(final) <- final[,1]
  final <- final[,-1]
  return(list(final,final_names))
}

pooler_check <- function(final,final_names, number_of_parts){
  final[final>1] <- 1
  if (number_of_parts == 1) {
    sums <- apply(final,2, sum)
    good_pools <- names(sums)[sums == nrow(final)]
    if (isEmpty(good_pools)) {
      return(FALSE)
    }
    else if (length(good_pools)>1) {
      message("Multiple suitable pools found. Starting scorring.")
      scores <- pooler_scoring(output_dir,good_pools)
      message("Calculated scores:")
      for(i in 1:length(good_pools)) message(paste(good_pools[i],": ",scores[i], sep="")) 
      return(good_pools)
    }
    else {
      message(paste("One suitable pool found: ",good_pools)); return(good_pools)
    }
  } else{
    all_comb <- combn(names(final),number_of_parts)
    control <- apply(all_comb,2, function(x, final) apply(final[,x],1,sum), final=final)
    best_pools_names <- all_comb[,apply(control,2, function(x) all(x==1))]
    combinations <- data.frame()
    if (!(isEmpty(best_pools_names))) {
      message("You can combine these pools: ")
      for(i in 1:length(best_pools_names)) message(best_pools_names[i]) 
      return(best_pools_names)
    } else {
      # primer assigned to both pools?
      control[control>1] <- 1
      best_pools_names <- all_comb[,apply(control,2, function(x) all(x==1))]
      
      if (!(isEmpty(best_pools_names))) {
        combinations <- combine_pools(final,final_names)
      } else return(FALSE)
    }
  }
}

pooler_scoring <- function(output_dir,names){
  system(paste('echo "Primer pooler scoring" > ', output_dir, "/pooler_output/poolfiles_score.txt",sep=""))
  for (i in 1:length(names)) {
    system(paste('echo "',names[i],'"', " >> " , output_dir, "/pooler_output/poolfiles_score.txt",sep=""))
    try(system(
        paste(primer_pooler, " ", output_dir, "/pooler_output/",names[i], " --print-bonds=1 ",
              ">> ", output_dir, "/pooler_output/poolfiles_score.txt", sep="")
      ))
      
  }
  
  scores <- readLines(paste(output_dir,"pooler_output/poolfiles_score.txt", sep="/"))
  where <- grep("poolfile\\d",scores)
  if (isEmpty(where)) where <- grep("Combination\\d+",scores)
  where <- c(where, length(scores))
  
  scores_sum <- c()
  for (i in 1:(length(where)-1)) {
    my_pfile <- scores[where[i]:where[i+1]]
    pool_score <- sum(as.integer(regmatches(my_pfile[grep("Score = ", my_pfile)],regexpr("\\d+",my_pfile[grep("Score = ", my_pfile)]))>2))
    scores_sum <- c(scores_sum, pool_score)
  }
  return(scores_sum)
}

pooler_picking <- function(final,final_names,success, number_of_parts,output_dir){
  final[final>1] <- 1
  if (number_of_parts == 1) {
    sums <- apply(final,2, sum)
    good_pools <- names(sums)[sums == nrow(final)]
    if (isEmpty(good_pools)) {
      message("No pool found"); return(NULL)
    }
    else if (length(good_pools)>1) {
      scores <- pooler_scoring(output_dir,good_pools)
      message("Calculated scores:")
      for(i in 1:length(good_pools)) message(paste(good_pools[i],": ",scores[i], sep="")) 
      return(good_pools)
    }
    else {
      message(paste("You can you this pool:",good_pools)); return(good_pools)
    }
  } else{
    # TO DO !!!
    all_comb <- combn(names(final),number_of_parts)
    control <- apply(all_comb,2, function(x, final) apply(final[,x],1,sum), final=final)
    best_pools_names <- all_comb[,apply(control,2, function(x) all(x==1))]
    
    if (!(isEmpty(best_pools_names))) return(best_pools_names) 
    else {
      # check for colision
      return(FALSE)
    }
  }
}

# Number of pools
primers <- names(readDNAStringSet(paste(output_dir,"primers.fasta", sep="/")))
number_of_parts <- regmatches(primers,regexpr("_\\d_p\\d-[RF]$",primers))
if (!(isEmpty(number_of_parts))){
  number_of_parts <- max(as.integer(gsub("_","",gsub("_p\\d-[FR]","",number_of_parts))))
} else number_of_parts <- 1

if (is.na(pools)) {
  primers <- regmatches(primers,regexpr("_p\\d-[RF]$",primers))
  pools <- max(as.integer(regmatches(primers,regexpr("\\d+",primers)))) + number_of_parts # +1 because of indexing and -1 because if we have two parts, we want to add 1
}

message(paste("Number of pools:",pools))

# create directory for pooler output
if (!(dir.exists((paste(output_dir, "pooler_output", sep = "/"))))) {
  dir.create(paste(output_dir, "pooler_output", sep="/"))
}

# run primer pooler 
success <- FALSE
attempts <- 1

while ((paste(success,collapse = " ") == FALSE) & (attempts <= 5)) {
  message(paste("Attempt No.", attempts))
  try(system(
    paste(primer_pooler ," --pools=", pools,",1,",output_dir,"/pooler_output/poolfile ", "--genome=fastGENdesigner_files/hg38.2bit ",output_dir, 
          "/primers.fasta 2> ", output_dir, "/pooler_output.txt", sep="")
  ))
  list_dataframes <- suppressWarnings(pooler_summary(output_dir))
  success <- pooler_check(list_dataframes[[1]],list_dataframes[[2]], number_of_parts)
  attempts <- attempts + 1
}

if (paste(success,collapse = " ") != FALSE) {
  excel_update(output_dir,list_dataframes[[1]],list_dataframes[[2]], success)
  message("Poolfiles and pooler output created")
  
} else {
  message("Giving up, trying to create some combinations")
}

#pooler_picking(list_dataframes[[1]],list_dataframes[[2]], success, number_of_parts, output_dir)
# message("Poolfiles and pooler output created")
# 
# 
# source("pooler_summary.R")
# 
# if(is.data.frame(which_to_test)) {
#   message("Scoring multiple combinations")
#   rownames(which_to_test) <- paste("Combination",rownames(which_to_test),sep="")
#   primers <- readDNAStringSet(paste(output_dir,"primers.fasta", sep="/"))
#   system(paste('echo "Primer pooler scoring" > ', output_dir, "/pooler_output/poolfiles_score.txt",sep=""))
#   for (i in 1:nrow(which_to_test)) {
#     temp_poolfile=tempfile()
#     writeXStringSet(primers[gsub("-[FR]","",names(primers)) %in% as.vector(t(which_to_test[i,]))],temp_poolfile)
#     system(paste('echo "Combination', i,' " >> ', output_dir, "/pooler_output/poolfiles_score.txt",sep=""))
#     
#     try(system(
#       paste(primer_pooler, " ", temp_poolfile, " --print-bonds=1 ",
#             ">> ", output_dir, "/pooler_output/poolfiles_score.txt", sep="")
#     ))
#     
#   }
#   source("pooler_scoring.R")
#   which_to_test["Score"] <- scores_sum
#   excel_update(output_dir = output_dir, combinations = which_to_test)
# } else message(paste("The best pools:",which_to_test))
# 
# #sink(file=NULL, type="output")
# #sink(file=NULL, type="message")
# 
# message("")
