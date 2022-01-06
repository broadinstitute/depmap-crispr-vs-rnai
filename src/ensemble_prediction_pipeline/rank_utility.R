

percentile_rank <- function(vec){
  
  vec[is.na(vec)] <- Inf
  vec <- rank(vec,ties.method = "max")
  vec <- vec / length(vec)
  
  return(vec)
}
