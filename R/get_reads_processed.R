

get_reads_processed.single <- function(metaphlan_profile.bfc){
  return(parse_number(readLines(metaphlan_profile.bfc, n = 3)[[3]]))
}

get_reads_processed.vec <- function(metaphlan_profiles.bfc){
  return(as.integer(sapply(metaphlan_profiles.bfc, function(x) readr::parse_number(readLines(x, n = 3)[[3]]))))
}
