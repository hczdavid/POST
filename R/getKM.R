#--------------------------------#
# Get kernel matrix functions
# Author:  David Huang
# Advisor: Jung-Ying Tzeng
# Date:    September 2020
#--------------------------------#

#---------------------------#
# Among OTU distance matrix
#---------------------------#

#' @title distanceMatrix
#' @param tree the phylogenetic tree
#' @param OTU the OTU table
#' @keywords internal
distanceMatrix <- function(tree,OTU) {
  Dmat       <- as.matrix(cophenetic(tree))
  sd_dMatrix <- sd(x = Dmat)
  otuname    <- colnames(OTU)
  Dmat       <- Dmat[otuname,otuname]
  return(list("dMatrix"    = Dmat,
              "sd_dMatrix" = sd_dMatrix))
}


#-----------------------#
# getkernel function
#-----------------------#
#' @title getkernel
#' @param dv distance matrix vector
#' @param OTU the OTU table
#' @param hx h value, c * sd
#' @keywords internal
getkernel <- function(dv, h, OTU){

  if(h == 0){
    Rmc <- exp(x = -dv^2/1e-10)
    }else{
    Rmc <- exp(x = -dv^2/h)
    }

  nsub <- nrow(OTU)
  Rmc  <- as.vector(Rmc)
  AD   <- matrix(0, nrow = nsub, ncol = nsub)
  for (ii in 2:nsub) {
    for (jj in 1:(ii-1)) {
      AD[jj,ii] = AD[ii,jj] = sqrt(sum((OTU[ii,]-OTU[jj,])^2*Rmc))
    }
  }

  ve     <- as.matrix(rep(1,nsub))
  mat    <- diag(nrow = nsub)-(ve%*%t(ve))/nsub
  kernel <- -1/2*mat%*%AD^2%*%mat
  return(list(kernel))
}


########################
#getkernellist function
########################
#' @title getkernellist
#' @param tree the phylogenetic tree
#' @param OTU the OTU table
#' @param cValues c vector
#' @keywords internal
getkernellist <- function(cValues, tree, OTU, verbose = TRUE) {

  OTU  <- as.matrix(OTU)
  nsub <- nrow(OTU)
  nOTU <- ncol(OTU)

  #add the pseudo-count
  OTU  <- OTU + 1e-6

  #CLR transformation
  for(i in 1:nsub){
    OTU[i,]=log(OTU[i,])-(1/nOTU)*sum(log(OTU[i,]))
  }

  dInfo      <- distanceMatrix(tree,OTU)
  kernellist <- list()
  for(w in 1:length(cValues)){
    if (verbose) cat("cValues", cValues[w], "\n")
    h <- cValues[w]*dInfo$sd_dMatrix
    kernellist[[w]] <- apply(dInfo$dMatrix,2,function(x){getkernel(x, h=h, OTU=OTU)})
  }
  return(kernellist)
}





