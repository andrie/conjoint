# 
# Author: Andrie
###############################################################################
#
#designMaxdiff <- function(choices, nAttributes, blocksize, 
#    versions=1, nrepeats=100, AttrNames=LETTERS(1:nAttributes)){
#  #choices <- 4
#  #nAttributes <- 14
#  #blocksize <- 7
#  #versions <- 2
#  bsize <- rep(rep(choices, blocksize), versions)
#  
#  tr   <- factor(1:nAttributes)
#  
#  designBlock
  

#' Create optimal block design for maxdiff scaling
#' 
#' This is a wrapper around \code{\link[AlgDesign]{optBlock}} to create partially balanced incomplete block designs.
#'
#' @param nAttributes The number of attributes to be tested
#' @param AttrNames A character vector with names for each attribute.
#' @param choices The number of choice tasks in each screen (typically about 4)
#' @param sets The desired number of sets (i.e. screens in market research terms)
#' @param versions The number of different versions in the design.  Each respondent will only see a single version, but the versions will be randomised across different respondents.
#' @param nrepeats The number of time the algorithm will repeat to find an optimal solution
#' @export 
designMaxdiff <- function(
    nAttributes, 
    AttrNames=LETTERS(1:nAttributes), 
    choices, 
    sets, 
    versions=1, 
    nrepeats=100
){
  
  tr   <- factor(1:nAttributes)
  designBlock(data=tr, choices=choices, sets=sets, 
      versions=versions, nrepeats=nrepeats)
}

###############################################################################

#' Create optimal block design for conjoint analysis scaling
#' 
#' This is a wrapper around \code{\link[AlgDesign]{optBlock}} to create partially balanced incomplete block designs.
#'
#' @param fractionalDesign The result of \code{\link{designConjoint}}
#' @param choices The number of choice tasks in each screen (typically about 4)
#' @param sets The desired number of sets (i.e. screens in market research terms)
#' @param versions The number of different versions in the design.  Each respondent will only see a single version, but the versions will be randomised across different respondents.
#' @param nrepeats The number of time the algorithm will repeat to find an optimal solution
#' @export 
designConjointBlock <- function(
    fractionalDesign, 
    choices, 
    sets, 
    versions=1, 
    nrepeats=100
){
  designBlock(data=fractionalDesign, choices=choices, sets=sets, 
      versions=versions, nrepeats=nrepeats)
}
  
designBlock <- function(
    data, 
    choices, 
    sets, 
    versions, 
    nrepeats){  
  bsize <- rep(rep(choices, sets), versions)
  crit <- ifelse(
      versions==1,
      "D",  # Best for single version?
      "Dp" # Best for multiple versions?
  )
  PBIB <- optBlock(~., withinData=data, blocksizes=bsize,
      nRepeats=nrepeats, criterion=crit)
  cat("\n\n=== Level by level table ===\n\n")
  print(as.table(crossprod(table(rep(1:(sets*versions), each=choices), PBIB$rows))),
      zero.print=".")
  
  for (i in 1:versions){
    cat("\n\n=== Set ",i," ===\n\n")
    print(as.table(crossprod(
            table(
                rep(1:sets, each=choices),
                PBIB$rows[(((i-1)*sets*choices)+1):(i*sets*choices)])
        )), zero.print='.')
  }
  cat("\n\n")
  
  PBIB$rows
  md <- t(array(PBIB$rows, dim=c(choices, sets*versions)))
  colnames(md) <- paste("Item", 1:choices, sep="")
  #ortho_md <- orthogonalBlock(md)
  ortho_md <- md
  ret <- data.frame(
      version=rep(1:versions, each=sets),
      set=rep(1:sets, versions),
      ortho_md
  )
  return(ret)
}

