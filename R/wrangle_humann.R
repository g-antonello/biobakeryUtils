wrangle_humann <- function(input_table, input_type = "ko"){
  
  if(input_type == "ko"){
  description.list <- strsplit(input_table[[1]], "\\ ") 
  
  # extract KO information
  ko_only <- gsub(":", "", sapply(description.list, "[", 1), fixed = TRUE)
  # extract EC information
  EC_only <- gsub("\\[|\\]", "", sapply(description.list, function(x) x[length(x)])) 
  # extract enzyme name information
  enz_name <- sapply(description.list, function(x) 
      tryCatch(x[2:(length(x) - 1)], 
               error = function(e) return(x))) |> 
    sapply(paste, collapse = " ")
  
  # build a data.frame with these data
  description.df <- data.frame(
    KO = ko_only,
    EC = EC_only,
    enzyme_name = enz_name
  )
  rownames(description.df) <- description.df[[1]]
  }
 
  if(input_type == "pathways"){
    description.list <- strsplit(input_table[[1]], "\\:\\ ") 
    
    description.df <- data.frame(
      pathway_short = sapply(description.list, "[", 1),
      pathway_long = sapply(description.list, "[", 2)
    ) 
  rownames(description.df) <- description.df[[1]]
   
  }
  
  # clean data matrix
  input_table_matrix <- as.matrix(input_table[,2:ncol(input_table)])
  rownames(input_table_matrix) <- rownames(description.df)
  
  
  # return results
  return(
    list(
      profiles = input_table_matrix,
      description = description.df
    )
  )
}