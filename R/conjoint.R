## Conjoint analysis experimental design
## 
## Author: Andrie
################################################################################

#' Create conjoint analysis experimental design.
#' 
#' @param nLevels A numeric vector specifying the number of levels for each attribute
#' @param extraTrials the number of additional trials to use in experiment
#' @param levelNames A character vector with Attribute names
#' @param ... Other arguments passed to \code{\link[AlgDesign]{optFederov}}
#' @export
#' @importFrom AlgDesign optFederov optBLock
#' @importFrom base LETTERS
#' @examples
#' designConjoint(c(6, 6, 3, 2))
designConjoint <- function(nLevels=c(4,3,2,1), extraTrials=3, levelNames=LETTERS[1:length(nLevels)], ...){
  #nLevels <- c(8, 6, 3, 3, 3, 2, 2, 2, 2, 2, 2) # vector with number of levels for each attribute
  #nLevels <- c(4, 4, 2, 2, 3, 2, 2, 2, 3, 2, 2, 3, 2, 2, 2, 2, 3, 2, 4)
  
  n <- length(levelNames)
  
  # Determine minimum number of trials to calculate first order effects
  minTrials <- sum(nLevels) + 1 - length(nLevels)
  
  cat("Number of attributes: ", length(nLevels),"\n")
  cat("Number of levels: ", sum(nLevels),"\n")
  cat("Minimum number of trials for 1st order effects: ", minTrials, "\n")
  cat("Number of stimuli for maximum efficiency: ", minTrials+extraTrials,"\n")
  
  # Generate the full factorial design, to use as input into Federov optimisation
  message("Generating full factorial design")
  factorialDesign <- gen.factorial(levels=nLevels, nVars=n, factors="all",
      varNames=levelNames)
  message("Estimating best fractional factorial design")
  fractionalDesign <- optFederov(~., data=factorialDesign, nRepeats=10, nTrials=minTrials + extraTrials, ...)
  
  cat("Design:\n\n")
  print(fractionalDesign$design)
  
  invisible(
      structure(
        list(
          nLevels          = nLevels,
          levelNames       = levelNames,
          factorialDesign  = factorialDesign,
          fractionalDesign = fractionalDesign$design
        ),
        class="conjoint"
    )
  )
}

#' Finds optimal conjoint design
#' 
#' @export 
bestDesignConjoint <- function(dat){
  stop("Function not yet defined")
  fractional <- function(extraTrials){
    optFederov(~., data=dat, nRepeats=10, nTrials=minTrials + extraTrials)
  } 
  
  effGe <- vapply(extraTrials, function(i)fractional(i)$Ge, 1)
  effGe
  #ggplot(data.frame(x=effD, y=effGe, trial=extraTrials), aes(x=x, y=y)) + geom_text(aes(label=trial))
  extraTrials <- extraTrials[which(effGe==max(effGe))]
  
  fractDesign <- fractional(extraTrials[1])
  
  cat("Efficiency: ", fractDesign$Ge,"\n",sep="")
  cat("Design:\n\n")
  print(fractDesign$design)
  return(null)
}


