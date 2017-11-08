
#########################################################
# Maxdiff scaling experimental design
# using partially balanced incomplete block design (PBIB)
#########################################################


#' Finds an orthogonal design for a given PBIB design.
#' 
#' @param x A partially balanced incomplete block (PBIB) experimental design
#' @param versions Number of versions of design
#' @param iter Number of iterations
#' @param restart TRUE or FALSE
#' @author Joris Meys
#' @export
#' @examples 
#' #Create test data
#' library(AlgDesign)
#' set.seed(12345)
#' choices <- 4
#' nAttributes <- 7
#' blocksize <- 7
#' bsize <- rep(choices, blocksize)
#' # Find optimal block design
#' PBIB <- optBlock(~., withinData=factor(1:nAttributes), blocksizes=bsize)
#' df <- data.frame(t(array(PBIB$rows, dim=c(choices, blocksize))))
#' colnames(df) <- paste("Item", 1:choices, sep="")
#' rownames(df) <- paste("Set", 1:nAttributes, sep="")
#' print(df)
#' set.seed(12345)
#' orthogonalBlock(df)
orthogonalBlock <- function(x , versions=nlevels(factor(x$version)), 
    iter=nrow(x)/length(unique(versions))+2, restart=TRUE){
  
  savedCols <- x[, c("version", "set")]
  savedNames <- names(x)
  x <- x[, -match(c("version", "set"), names(x))]
  

  balanceMatrix <- function(x){
    t(sapply(sort(unique(unlist(x))), function(i)colSums(x==i)))
  }
  
  balanceScore <- function(x){
    sum((1-x)^2)
  }
  
  opt_block <- function(x, iter, restart){
    # transform rows to list
    sets <- unlist(apply(x, 1, list), recursive=FALSE)
    nsets <- NROW(x)
    # C0 contains all possible design points
    C0 <- lapply(sets, combinat::permn)
    n <- gamma(NCOL(x)+1)
    
    # starting point
    id <- sample(1:n, nsets)
    Sol <- sapply(1:nsets, function(i)C0[[i]][id[i]])
    
    IT <- iter
    score_best <- Inf
    out_best <- x
    # other iterations
    while(IT > 0){
      message(paste(rep(".", IT), collapse=""))
      for(i in 1:nsets){
        nn <- 1:n
        scores <- sapply(nn, function(p){
              tmp <- Sol
              tmp[[i]] <- C0[[i]][[p]]
              w <- balanceMatrix(do.call(rbind, tmp))
              balanceScore(w)
            })
        idnew <- nn[which.min(scores)]
        Sol[[i]] <- C0[[i]][[idnew]]
        
      }
      #Check if score is 0
      out <- as.data.frame(do.call(rbind, Sol))
      score <- balanceScore(balanceMatrix(out))
      if(score < score_best){
        out_best <- out
        score_best <- score
      }  
      if (score==0) {break}
      IT <- IT - 1
      
      # If asked, restart
      if(IT==0 & restart){
        restart <- FALSE
        message("restarting")
        id <- sample(1:n, nsets)
        Sol <- sapply(1:nsets, function(i)C0[[i]][id[i]])
        IT <- iter
      }
    }
    message(paste("Orthogonal block: score = ", score_best))
    out_best
  }
  
  stopifnot(require(combinat))
  
  if(identical(versions,1)){
    ret <- opt_block(x, iter, restart)
  } else {
    ret <- lapply(unique(versions), function(xt)opt_block(x[xt==versions, ], iter=iter, restart=restart))
    ret <- do.call(rbind, ret)
  }
  
  ret <- data.frame(savedCols, ret)
  names(ret) <- savedNames
  ret
  
}





#orthogonalBlock <- function(x, iter=10, restart=FALSE){
#  iter_orig <- iter
#  
#  balanceMatrix <- function(x){
#    t(sapply(sort(unique(unlist(x))), function(i)colSums(x==i)))
#  }
#  
#  balanceScore <- function(x){
#    sum((1-x)^2)
#  }
#  
#  stopifnot(require(combinat))
#  # transform rows to list
#  sets <- unlist(apply(x, 1, list), recursive=FALSE)
#  nsets <- NROW(x)
#  # C0 contains all possible design points
#  C0 <- lapply(sets, permn)
#  n <- gamma(NCOL(x)+1)
#  
#  # starting point
#  id <- sample(1:n, nsets)
#  Sol <- sapply(1:nsets, function(i)C0[[i]][id[i]])
#
#  score_best <- Inf
#  out_best <- x
#  
#  # other iterations
#  while(iter > 0){
#    message(iter)
#    for(i in 1:nsets){
#      nn <- (1:n)
#      scores <- sapply(nn,function(p){
#            tmp <- Sol
#            tmp[[i]] <- C0[[i]][[p]]
#            w <- balanceMatrix(do.call(rbind, tmp))
#            balanceScore(w)
#          })
#      id[i] <- nn[which.min(scores)]
#      Sol[[i]] <- C0[[i]][[id[i]]]
#      
#    }
#    #Check if score is 0
#    out <- as.data.frame(do.call(rbind, Sol))
#    score <- balanceScore(balanceMatrix(out))
#    
#    if(score < score_best) out_best <- out
#    message(score)
#    if (score==0) {break}
#    iter <- iter - 1
#    
#    # If asked, restart
#    if(iter==0 && restart){
#      restart <- FALSE
#      id <- sample(1:n, nsets)
#      Sol <- sapply(1:nsets, function(i)C0[[i]][id[i]])
#      iter <- iter_orig
#    }
#  }
#  out_best
#}








