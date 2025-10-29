wrangle_humann <- function(input_table, input_type = "ko"){
  
  if(input_type == "ko"){
  featureName <- input_table[[1]]
  
  # extract KO information
  ko_only <- sapply(strsplit(featureName, ": ", fixed = TRUE), "[", 1)
  # extract EC information
  EC_only <- tryCatch(sapply(strsplit(featureName, "[", fixed = TRUE), function(x) ifelse(length(x) > 1, gsub("]", "", x[length(x)]), NA)), error = function(e) return(NA))
  # extract enzyme name information
  enz_name <- tryCatch(sapply(strsplit("K00383: glutathione reductase (NADPH) [EC:1.8.1.7]", "\\:\\ |\\ \\["), "[", 2), error = function(e) return(NA))
  
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
      pathway = sapply(description.list, "[", 1),
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