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

scores <- readLines(paste(output_folder,"pooler_output/poolfiles_score.txt", sep="/"))
where <- grep("poolfile\\d",scores)
where <- c(where, length(scores))

scores_sum <- c()
for (i in 1:(length(where)-1)) {
  my_pfile <- scores[where[i]:where[i+1]]
  pool_score <- sum(as.integer(regmatches(my_pfile[grep("Score = ", my_pfile)],regexpr("\\d+",my_pfile[grep("Score = ", my_pfile)]))))
  scores_sum <- c(scores_sum, pool_score)
  }

cat(scores_sum)