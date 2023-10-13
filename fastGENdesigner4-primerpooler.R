# STEP4 -PRIMER POOLER
cat("Starting primer pooler\n")

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

excel_update <- function(output_dir,final=NULL,final_names=NULL, best_pools = NULL, combinations=NULL,best_combination_list,output_setting,input_file, final_bed) {
  files <- list.files(path = paste(output_dir,"pooler_output", sep="/"), pattern="poolfile\\d", full.names = T) 
  # adding dataframes into Excel file
  wb <- loadWorkbook(paste(output_dir,"fastGENdesigner-output.xlsx", sep ="/"))
  
  # adjusting primer properties - excluded primer pairs are marked as red
  before <- read.xlsx(wb, "Primer properties",colNames=TRUE)$Name
  after <- names(readDNAStringSet(paste(output_dir,"primers.fasta", sep ="/")))
  writeData(wb,"Primer properties",data.frame(Selected=as.integer(gsub("-[RF]","",before) %in% final_bed)),startCol = 10)
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
                                       sheet_name = "All attempts") 
    writeData(wb, sheet = sheet_name, combinations,rowNames=TRUE,colNames = TRUE)
    #writeData(wb, sheet = sheet_name, best_combination,rowNames=TRUE, startRow =nrow(combinations) + 4)
    setColWidths(wb,  sheet_name, cols = 1:(ncol(combinations)+1), widths = "auto")
    sheet_name <- add_unique_worksheet(wb, file = paste(output_dir,"fastGENdesigner-output.xlsx", sep ="/"),
                                       sheet_name = "Combinations") 
    writeData(wb, sheet = sheet_name, best_combination_list[[1]], rowNames = TRUE, colNames = TRUE)
    setColWidths(wb,  sheet_name, cols = 1:3, widths = "auto")
    addStyle(wb, sheet_name,cols=1:3, rows=which(rownames(best_combination_list[[1]]) %in% unlist(best_combination_list[[2]]))+1,
             style = createStyle(bgFill= "green"), gridExpand=TRUE, stack = FALSE)
    
  } else{
    cat("Primers distribution saved to Excel file: fastGENdesigner-output.xlsx\n")
  }
  sheet_name <- add_unique_worksheet(wb,file = paste(output_dir,"fastGENdesigner-output.xlsx",sep="/"), sheet_name = "Settings")
  writeData(wb,sheet=sheet_name,output_setting,rowNames=TRUE,colNames = FALSE)
  writeData(wb,sheet=sheet_name,input_file,startCol=ncol(output_setting) + 3)
  setColWidths(wb,  sheet_name, cols = 1:2, widths = "auto")
  
  saveWorkbook(wb,paste(output_dir,"fastGENdesigner-output.xlsx", sep ="/"), overwrite = TRUE)
  
}

pooler_combinations <- function(final, final_names, best_pools_names){
  LETTERS702 <- sapply(LETTERS, function(x) paste0(x, LETTERS))
  my_LETTERS <- c(LETTERS702,sapply(LETTERS702, function(x) paste0(x, LETTERS)))
  best_pools_names <- as.data.frame(best_pools_names)
  combinations <- data.frame()
  #final_names$primer_pairs <- gsub("_p$","", final_names$primer_pairs)
  final <- final[unique(unlist(best_pools_names))]
  rownames(final) <- substring(rownames(final),1,nchar(rownames(final))-2)
  final_names <- final_names[unique(unlist(best_pools_names))]
  rownames(final_names) <- rownames(final)
  control <- apply(final,1,sum)
  numbers_of_combs <- c() 
  for (j in 1:ncol(best_pools_names)){
    for (pool_name in best_pools_names[,j]){
      my_list <- list()
      for (name in names(control)) {
        pridat <- unlist(strsplit(final_names[name,pool_name], split=" "))
        if (paste(pridat, collapse = "") != "-") pridat <- paste(name,pridat,sep="_")
        if (control[name] > 1) pridat <- c(pridat,"-")
        my_list[[name]] <-  pridat
      }
      comb <- expand.grid(my_list)
      # vymazat nezmyselne kombinacie - musi obsahovat polovicu +- 3 pary
      balance_ind <- apply(comb,1,function(x) sum(x == "-"))
      #comb <- comb[(!((balance_ind <= round((nrow(final)/2)) -2) | (balance_ind >= round((nrow(final)/2)) +2))),]
      numbers_of_combs <- c(numbers_of_combs,nrow(comb))
    }}
  best_pools_names <- best_pools_names[,which.min(rowSums(matrix(numbers_of_combs, ncol = 2, byrow = TRUE)))]
  for (pool_name in best_pools_names){
    my_list <- list()
    for (name in names(control)) {
      pridat <- unlist(strsplit(final_names[name,pool_name], split=" "))
      if (paste(pridat, collapse = "") != "-") pridat <- paste(name,pridat,sep="_")
      if (control[name] > 1) pridat <- c(pridat,"-")
      my_list[[name]] <-  pridat
    }
    comb <- expand.grid(my_list)
    # vymazat nezmyselne kombinacie - musi obsahovat polovicu +- 1
    balance_ind <- apply(comb,1,function(x) sum(x == "-"))
    comb <- comb[(!((balance_ind <= round((nrow(final)/2)) -2) | (balance_ind >= round((nrow(final)/2)) +2))),]
    if (nrow(comb)>1) {
      my_names <- c()
      for (i in 1:nrow(comb)) my_names <- c(my_names,paste(pool_name,my_LETTERS[i]))
    } 
    else my_names <- name
    rownames(comb) <- my_names
    combinations <- rbind(combinations,comb)
  }
  return(combinations)
}

  #if (!((ncol(best_pools_names)>1) | (TRUE %in% (final[unique(unlist(best_pools_names))] > 1)))) return(NULL)
  # for (pool in unique(unlist(best_pools_names))) {
  #     my_list <- list()
  #     for (seq in final_names$primer_pairs) {
  #       pridat <- unlist(strsplit(final_names[seq,pool], split=" "))
  #       if (paste(pridat, collapse = "") != "-") pridat <- paste(seq,pridat,sep="_")
  #       my_list[[seq]] <-  pridat
  #     }
  #     comb <- expand.grid(my_list)
  #     if (nrow(comb)>1) {
  #       my_names <- c()
  #       for (i in 1:nrow(comb)) my_names <- c(my_names,paste(pool,LETTERS[i]))
  #     }
  #     else my_names <- pool
  #     rownames(comb) <- my_names
  #     combinations <- rbind(combinations,comb)
  # }
#  return(combinations)
#}

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
  best_pools_names <- as.data.frame(best_pools_names)
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
  #best_pools_names <- as.data.frame(best_pools_names)
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

pooler_check <- function(final,final_names, number_of_parts,attempt){
  final[final>1] <- 1
  if (number_of_parts == 1) {
    sums <- apply(final,2, sum)
    good_pools <- names(sums)[sums == nrow(final)]
    if (isEmpty(good_pools)) {
      return(FALSE)
    }
    else if (length(good_pools)>1) {
      cat("Multiple suitable pools found. Starting scorring.\n")
      scores <- pooler_scoring(output_dir,good_pools)
      cat("Calculated scores:\n")
      for(i in 1:length(good_pools)) cat(paste(good_pools[i],": ",scores[i], "\n", sep="")) 
      return(good_pools)
    }
    else {
      cat(paste("One suitable pool found: ",good_pools, "\n")); return(good_pools)
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
      cat("You can combine these pools: \n")
      for(i in 1:length(best_pools_names)) cat(paste(best_pools_names[i], "\n",)) 
      return(best_pools_names)
    } else {
      if (attempt <= 5) return(FALSE) 
      # SOME COMBINATIONS CAN BE MADE
      control[control>1] <- 1
      best_pools_names <- as.data.frame(all_comb[,apply(control,2, function(x) all(x==1))])
      # OVERLAP CHECK
      best_pools_names <- overlap_check(final,final_names,best_pools_names)
      if (isEmpty(best_pools_names)) return(FALSE)
      cat("You can combine these pools: \n")
      for(i in 1:length(best_pools_names)) cat(paste(best_pools_names[i],"\n")) 
      return(best_pools_names)
      #return(FALSE)
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
    all_scores <- as.integer(regmatches(my_pfile[grep("Score = ", my_pfile)],regexpr("\\d+",my_pfile[grep("Score = ", my_pfile)])))
    pool_score <- sum(all_scores[all_scores>2])
    scores_sum <- c(scores_sum, pool_score)
  }
  
  scores <- readLines(paste(output_dir,"pooler_output/poolfiles_dG.txt", sep="/"))
  where <- grep("poolfile\\d",scores)
  where <- c(where, length(scores))
  
  scores_sum_dG <- c()
  for (i in 1:(length(where)-1)) {
    my_pfile <- scores[where[i]:where[i+1]]
    all_scores <- as.numeric(gsub(",",".",regmatches(my_pfile[grep("dG = ", my_pfile)],regexpr("-\\d+,*\\d+",my_pfile[grep("dG = ", my_pfile)]))))
    #pool_score <- sum(all_scores[all_scores< -2])
    pool_score <- sum(all_scores)
    scores_sum_dG <- c(scores_sum_dG, pool_score)
  }
  
  return(list(scores_sum,scores_sum_dG))
}

pooler_scoring <- function(output_dir,names=NULL){
  system(paste('echo "Primer pooler scoring" > ', output_dir, "/pooler_output/poolfiles_score.txt",sep=""))
  if (is.null(names)) {
    names <- list.files(paste(output_dir,"pooler_output", sep="/"), pattern="poolfile\\d+.txt")
  }
  temp_file <- tempfile()
  for (i in 1:length(names)) {
    # SCORE
    system(paste('echo "',names[i],'"', " >> " , output_dir, "/pooler_output/poolfiles_score.txt",sep=""))
    try(system(
        paste(primer_pooler, " ", output_dir, "/pooler_output/",names[i], " --print-bonds=1 ",
              ">> ", output_dir, "/pooler_output/poolfiles_score.txt 2> ", temp_file, sep="")
      ))
    
    # dG
    system(paste('echo "',names[i],'"', " >> " , output_dir, "/pooler_output/poolfiles_dG.txt",sep=""))
    try(system(
      paste(primer_pooler, " ", output_dir, "/pooler_output/",names[i], 
            " --print-bonds=-2 --dg=333,2,1,1 ",">> ", output_dir, "/pooler_output/poolfiles_dG.txt 2> ",temp_file, sep="")
    ))
  }
  unlink(temp_file)
  scores_sum <- scores_summing(output_dir)
  names(scores_sum) <- c("Score","dG")
  return(scores_sum) 
  
}

combinations_scoring <- function(output_dir, combinations){
  cat("Scoring multiple combinations - this can take a while .. \n")
  primers <- readDNAStringSet(paste(output_dir,"primers.fasta", sep="/"))
  system(paste('echo "Primer pooler scoring" > ', output_dir, "/pooler_output/poolfiles_score.txt",sep=""))
  
  temp_file <- tempfile()
  for (i in 1:nrow(combinations)) {
    temp_poolfile=tempfile()
    writeXStringSet(primers[gsub("-[FR]","",names(primers)) %in% as.vector(t(combinations[i,]))],temp_poolfile)
    
    #SCORE
    system(paste('echo "', row.names(combinations)[i],'" >> ', output_dir, "/pooler_output/poolfiles_score.txt",sep=""))
    try(system(
      paste(primer_pooler, " ", temp_poolfile, " --print-bonds=1 ",
            ">> ", output_dir, "/pooler_output/poolfiles_score.txt 2> ", temp_file, sep="")
    ))
    
    system(paste('echo "',row.names(combinations)[i],'"', " >> " , output_dir, "/pooler_output/poolfiles_dG.txt",sep=""))
    try(system(
      paste(primer_pooler, " ", temp_poolfile, " --print-bonds=-2 --dg=333,2,1,1 ",
            ">> ", output_dir, "/pooler_output/poolfiles_dG.txt 2> ",temp_file, sep="")
    ))
    
  }
  unlink(temp_poolfile)
  scores <- scores_summing(output_dir)
  combinations["Score"] <- scores[[1]]
  combinations["dG"] <- scores[[2]]
  return(combinations)
}

combinations_choosing <- function(combinations,best_pools_names) {
  final_score <- c()
  final_score_dG <- c()
  my_names <- c()
  
  #for (i in 1:ncol(best_pools_names)) {
  #my_list <- list()
  #for (j in 1:length(best_pools_names)) {
  #  my_pool <- best_pools_names[j]
  #  my_list[[my_pool]] <- grep(my_pool,rownames(combinations), value=TRUE)
  #}
  #combs <- expand.grid(my_list)
  first_pool_comb <- combinations[grepl(best_pools_names[1],rownames(combinations)),]
  second_pool_comb <- combinations[!(grepl(best_pools_names[1],rownames(combinations))),]
  n_col <- ncol(combinations)-2 
  combs <- data.frame()
  for (r in 1:nrow(first_pool_comb)){
    a <- gsub("_p\\d+$","",second_pool_comb[,1]) != gsub("_p\\d+$","",first_pool_comb[r,1])
    for (c in 2:n_col){
      a <- (a == TRUE) & (gsub("_p\\d+$","",second_pool_comb[,c]) != gsub("_p\\d+$","",first_pool_comb[r,c]))
    }
    if (nrow(second_pool_comb[a,]) == 0) next
    combs <- rbind(combs,t(data.frame(c(rownames(first_pool_comb[r,]),rownames(second_pool_comb[a,])))))
  }
  rownames(combs) <- 1:nrow(combs)
  
  # keep only those combinations that are valid 
  #delete_ind <- c()
  
  #for (c in 1:nrow(combs)){
  #  control <- apply(combinations[unlist(combs[c,]),1:n_col],2, function(x) sum(grepl("-",x)))
  #  if (0 %in% control | 2 %in% control) delete_ind <- c(delete_ind,c)
  #}
  #  combs <- combs[-delete_ind,]
  for (c in 1:nrow(combs)) {
    final_score <- c(final_score,sum(combinations[unlist(combs[c,]),"Score"]))
    final_score_dG <- c(final_score_dG,sum(combinations[unlist(combs[c,]),"dG"]))
    my_names <- c(my_names,paste(unlist(combs[c,]), collapse = " + "))
  }
  final_scores <- data.frame(Score = final_score, dG = final_score_dG,row.names = my_names)
  # finding MINIMUM
  # first dG , second Score
  
  best_combination <- final_scores[final_scores[,"dG"] == max(final_scores[,"dG"]),]
  best_combination <- rownames(best_combination[best_combination[,"Score"] == min(best_combination[,"Score"]),])
  cat(paste("The best combination is: ",unlist(best_combination),". You can find it in the OUTPUT EXCEL TABLE.",  "\n",sep = ""))
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

cat(paste("Number of pools:",pools, "\n"))

# create directory for pooler output
if (!(dir.exists((paste(output_dir, "pooler_output", sep = "/"))))) {
  dir.create(paste(output_dir, "pooler_output", sep="/"))
}

# run primer pooler 
success <- FALSE
attempts <- 1

while ((paste(success,collapse = " ") == FALSE) & (attempts <= 5)) {
  cat(paste("Attempt No.", attempts, "\n"))
  try(system(
    paste(primer_pooler ," --pools=", pools,",1,",output_dir,"/pooler_output/poolfile ", "--genome=fastGENdesigner_files/hg38.2bit ",output_dir, 
          "/primers.fasta 2> ", output_dir, "/pooler_output/pooler_output.txt", sep="")
  ))
  file.remove(grep("overlap-report.\\d+.txt",list.files(), value=TRUE))
  list_dataframes <- suppressWarnings(pooler_summary(output_dir))
  success <- pooler_check(list_dataframes[[1]],list_dataframes[[2]], number_of_parts,attempts)
  attempts <- attempts +1
}

# USPECH PRE 1 POOL
if (paste(success,collapse = " ") != FALSE) {
  cat("Poolfiles and pooler output created\n")
  if (length(unlist(success)) == 1) combinations <- NULL
  else combinations <- pooler_combinations(list_dataframes[[1]],list_dataframes[[2]], success) # <---- problem
  best_combination_list <- NULL
  if (is.data.frame(combinations)) {
    combinations <- combinations_scoring(output_dir,combinations)
    best_combination_list <- combinations_choosing(combinations,success)
  }
  
  #excel_update(output_dir,list_dataframes[[1]],list_dataframes[[2]], success,combinations, best_combination_list[[1]],output_setting)
} else {
  # SKUSIT ZVYSIT POCET POOLOV
  s <- "s"
  if (number_of_parts == 1) s <- ""
  cat(paste("Cannot create ", number_of_parts, " pool",s,". Incrising number of pools.", "\n", sep = ""))
  attempts <- 6
  number_of_parts <- number_of_parts+1
  success <- pooler_check(list_dataframes[[1]],list_dataframes[[2]], number_of_parts,attempts)
  
  if (paste(success,collapse = " ") != FALSE) {
    combinations <- pooler_combinations(list_dataframes[[1]],list_dataframes[[2]], success)
    print(nrow(combinations))
    if (nrow(combinations)>2000) success = FALSE
  }
  
  while ((paste(success,collapse = " ") == FALSE) & (attempts <= 15)) {
    cat(paste("Attempt No.", attempts, "\n"))
    try(system(
      paste(primer_pooler ," --pools=", pools,",1,",output_dir,"/pooler_output/poolfile ", "--genome=fastGENdesigner_files/hg38.2bit ",output_dir, 
            "/primers.fasta 2> ", output_dir, "/pooler_output/pooler_output.txt", sep="")
    ))
    file.remove(grep("overlap-report.\\d+.txt",list.files(), value=TRUE))
    list_dataframes <- suppressWarnings(pooler_summary(output_dir))
    success <- pooler_check(list_dataframes[[1]],list_dataframes[[2]], number_of_parts,attempts)
    combinations <- pooler_combinations(list_dataframes[[1]],list_dataframes[[2]], success)
    print(nrow(combinations))
    if (nrow(combinations)>2500) success = FALSE
    else success <- unique(regmatches(rownames(combinations), regexpr("poolfile\\d+.txt", rownames(combinations))))
    attempts <- attempts +1
  }
    cat("Poolfiles and pooler output created\n")
    best_combination_list <- NULL
    if (is.data.frame(combinations)) {
      combinations <- combinations_scoring(output_dir,combinations)
      best_combination_list <- combinations_choosing(combinations,success)
    }
}

# SKOROVAT AJ ORIGINALNE POOLY
scores <- pooler_scoring(output_dir)
scores <- t(data.frame(scores))
colnames(scores) <- colnames(list_dataframes[[1]])
list_dataframes[[1]] <- rbind(list_dataframes[[1]], scores)
list_dataframes[[2]] <- rbind(list_dataframes[[2]], scores)
