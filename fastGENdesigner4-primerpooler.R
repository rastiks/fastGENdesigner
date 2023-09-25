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

excel_update <- function(output_dir,final=NULL,final_names=NULL, best_pools = NULL, combinations=NULL,best_combination) {
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
    
    addStyle(wb, sheet_name,rows=1:(nrow(final)+1), cols=which(colnames(final) %in% unlist(best_pools)) + 1,
             style = createStyle(bgFill= "green"), gridExpand=TRUE, stack = FALSE)
    
    sheet_name <- add_unique_worksheet(wb, file = paste(output_dir,"fastGENdesigner-output.xlsx", sep ="/"), 
                                       sheet_name = paste("PrimerDistributionNames-",length(files), sep =""))
    writeData(wb, sheet = sheet_name, final_names)
    setColWidths(wb,  sheet_name, cols = 1, widths = "auto")
    addStyle(wb, sheet_name,rows=1:(nrow(final_names)+1), cols=which(colnames(final_names) %in% unlist(best_pools)) ,
             style = createStyle(bgFill= "green"), gridExpand=TRUE, stack = FALSE)
  }
  
  if (!(is.null(combinations)) & is.data.frame(combinations)) {
    sheet_name <- add_unique_worksheet(wb, file = paste(output_dir,"fastGENdesigner-output.xlsx", sep ="/"),
                                       sheet_name = "Combinations") 
    writeData(wb, sheet = sheet_name, combinations,rowNames=TRUE)
    writeData(wb, sheet = sheet_name, best_combination_list,rowNames=TRUE, startRow =nrow(combinations) + 3)
    setColWidths(wb,  sheet_name, cols = 1:(ncol(combinations)+1), widths = "auto")
  } else{
    message("Primers distribution saved to Excel file: fastGENdesigner-output.xlsx")
  }
  
  saveWorkbook(wb,paste(output_dir,"fastGENdesigner-output.xlsx", sep ="/"), overwrite = TRUE)
  
}

pooler_combinations <- function(final, final_names, best_pools_names){
  best_pools_names <- as.data.frame(best_pools_names)
  combinations <- data.frame()
  final_names$primer_pairs <- gsub("_p$","", final_names$primer_pairs)
  rownames(final_names) <- final_names$primer_pairs
  if (!((ncol(best_pools_names)>1) | (TRUE %in% (final[unique(unlist(best_pools_names))] > 1)))) return(NULL)
  for (pool in unique(unlist(best_pools_names))) {
      my_list <- list()
      for (seq in final_names$primer_pairs) {
        pridat <- unlist(strsplit(final_names[seq,pool], split=" "))
        if (paste(pridat, collapse = "") != "-") pridat <- paste(seq,pridat,sep="_")
        my_list[[seq]] <-  pridat
      }
      comb <- expand.grid(my_list)
      if (nrow(comb)>1) {
        my_names <- c()
        for (i in 1:nrow(comb)) my_names <- c(my_names,paste(pool,LETTERS[i]))
      }
      else my_names <- pool
      rownames(comb) <- my_names
      combinations <- rbind(combinations,comb)
  }
  return(combinations)
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

combine_pools <- function(final,final_names,best_pools_names) {
  for (i in 1:ncol(best_pools_names)) {
    groups <- gsub("_\\d_p$","",final_names$primer_pairs)
    final_names$primer_pairs <- groups
    groups <- unique(groups)
    if (TRUE %in% ((apply(final[final_names$primer_pairs==groups[i],],2, sum))>1)) {
      } # check if there is an overlap
    print(best_pools_names[,i])
  }
  
}

overlap_check <- function(final,final_names, best_pools_names){
  best_pools_names <- as.data.frame(best_pools_names)
  groups <- gsub("_\\d_p$","",final_names$primer_pairs)
  final_names$primer_pairs <- groups
  groups <- unique(groups)
  overlap <- c()
  for (i in 1:ncol(best_pools_names)){
    best_df <- final[best_pools_names[,i]]
    best_df_names <- final_names[c("primer_pairs",best_pools_names[,i])]
    for (j in 1:length(groups)) {
      if ((TRUE %in% ((apply(best_df[best_df_names$primer_pairs==groups[j],],2, sum))>1))) {
        overlap <- c(overlap,i)}
    }
  }
  if (!is.null(overlap)) best_pools_names <- best_pools_names[,-overlap]
  return(best_pools_names)
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
    best_pools_names <- as.data.frame(all_comb[,apply(control,2, function(x) all(x==1))])
    combinations <- data.frame()
    if (!(isEmpty(best_pools_names))) {
      # OVERLAP CHECK
      best_pools_names <- overlap_check(final,final_names,best_pools_names)
      if (isEmpty(best_pools_names)) return(FALSE)
      message("You can combine these pools: ")
      for(i in 1:length(best_pools_names)) message(best_pools_names[i]) 
      return(best_pools_names)
    } else {
      # SOME COMBINATIONS CAN BE MADE
      control[control>1] <- 1
      best_pools_names <- all_comb[,apply(control,2, function(x) all(x==1))]
      # OVERLAP CHECK
      best_pools_names <- overlap_check(final,final_names,best_pools_names)
      if (isEmpty(best_pools_names)) return(FALSE)
      return(best_pools_names)
    }
  }
}

scores_summing <- function(output_dir) {
  scores <- readLines(paste(output_dir,"pooler_output/poolfiles_score.txt", sep="/"))
  where <- grep("poolfile\\d",scores)
  where <- c(where, length(scores))
  
  scores_sum <- c()
  for (i in 1:(length(where)-1)) {
    my_pfile <- scores[where[i]:where[i+1]]
    pool_score <- sum(as.integer(regmatches(my_pfile[grep("Score = ", my_pfile)],regexpr("\\d+",my_pfile[grep("Score = ", my_pfile)]))>2))
    scores_sum <- c(scores_sum, pool_score)
  }
  return(scores_sum)
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
  scores_sum <- scores_summing(output_dir)
  
}

combinations_scoring <- function(output_dir, combinations){
  message("Scoring multiple combinations")
  primers <- readDNAStringSet(paste(output_dir,"primers.fasta", sep="/"))
  system(paste('echo "Primer pooler scoring" > ', output_dir, "/pooler_output/poolfiles_score.txt",sep=""))
  
  for (i in 1:nrow(combinations)) {
    temp_poolfile=tempfile()
    writeXStringSet(primers[gsub("-[FR]","",names(primers)) %in% as.vector(t(combinations[i,]))],temp_poolfile)
    system(paste('echo "', row.names(combinations[i]),'" >> ', output_dir, "/pooler_output/poolfiles_score.txt",sep=""))
    
    try(system(
      paste(primer_pooler, " ", temp_poolfile, " --print-bonds=1 ",
            ">> ", output_dir, "/pooler_output/poolfiles_score.txt", sep="")
    ))
  }
  
  scores_sum <- scores_summing(output_dir)
  combinations["Score"] <- scores_sum
  return(combinations)
}

combinations_choosing <- function(combinations,best_pools_names) {
  final_score <- c()
  my_names <- c()
  
  for (i in 1:ncol(best_pools_names)) {
    my_list <- list()
    for (j in 1:nrow(best_pools_names)) {
      my_pool <- best_pools_names[j,i]
      my_list[[my_pool]] <- grep(my_pool,rownames(combinations), value=TRUE)
    }
    combs <- expand.grid(my_list)
    for (c in 1:nrow(combs)) {
      final_score <- c(final_score,sum(combinations[unlist(combs[c,]),"Score"]))
      my_names <- c(my_names,paste(unlist(combs[c,]), collapse = " + "))
    }
  }
  final_scores <- data.frame(final_score,row.names = my_names)
  best_combination <- rownames(final_scores)[which.min(unlist(final_scores))]
  return(list(final_scores,best_combination))
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
  message("Poolfiles and pooler output created")
  combinations <- pooler_combinations(list_dataframes[[1]],list_dataframes[[2]], success)
  if (is.data.frame(combinations)) {
    combinations <- combinations_scoring(output_dir,combinations)
    best_combination_list <- combinations_choosing(combinations,success)
  }
  
  excel_update(output_dir,list_dataframes[[1]],list_dataframes[[2]], success,combinations, best_combination_list[[1]])
  
} else {
  message("Giving up, trying to create some combinations")
}

