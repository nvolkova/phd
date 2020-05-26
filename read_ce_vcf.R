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

genoToCoverage <- function(vcf){
  rowSums(sapply(grep("F|R.Z",names(geno(vcf)), value=TRUE), function(field){
    geno(vcf)[[field]]
  }, simplify="array"), dims=2) # number of all reads covering a particular position (substitution)
}

altCounts <- function(vcf){
  altAllele <- as.character(unlist(alt(vcf)))
  alleleCounts <- sapply(grep("F|R.Z",names(geno(vcf)), value=TRUE), function(field){
    geno(vcf)[[field]]
  }, simplify="array")
  t(sapply(seq_along(altAllele), function(i){
    alleleCounts[i,"TUMOUR",paste0(c("F","R"),altAllele[i], "Z")]
  }))
}

classifyRearrangements <- function(data){
  x <- as.character(data[,1]) != as.character(data[,5]) # diff./same chromosome
  y <- data[,2] == "+" # 1 orientation
  z <- data[,6] == "-" # 2 orientation
  cls <- pmin(4 * x + 2 * y + z, 4) 
  cls <- factor(cls, levels=0:4, labels = c("HHI","TD","DEL","TTI","INT"))
  cls
  # head-to-head insertion
  # tail-to-tail insertion
  # TD - tandem duplication: 
  # DEL - deletion?
  # interchromosomal
}

classifyRearrangements_old <- function(data){
  x <- as.character(data[,1]) != as.character(data[,5]) # diff./same chromosome
  y <- data[,2] == "+" # 1 orientation
  z <- data[,6] == "+" # 2 orientation
  cls <- pmin(4 * x + 2 * y + z, 4) 
  cls <- factor(cls, levels=0:4, labels = c("HHI","TD","DEL","TTI","INT"))
  cls
  # head-to-head insertion
  # tail-to-tail insertion
  # TD - tandem duplication: 
  # DEL - deletion?
  # interchromosomal
}


read_ce_table <- function(file, ...) {
  out <- tryCatch(
    {
      read.delim(file, ...) 
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