# tiny supplementary functions

top <- function(x,n) {
  return(sort(x,decreasing=T)[1:n])
}

cosine <- function(x,y) {
  sum(x * y) / sqrt(sum(x**2)) / sqrt(sum(y**2))
}

read_ce_vcf <- function(file) {
  out <- tryCatch(
    {
      readVcf(file, genome="WBcel235") 
    },
    error=function(cond) {
      message(paste("\nBad file in: ", file))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("\nFile caused a warning: ", file))
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    }
  )    
  return(out)
}