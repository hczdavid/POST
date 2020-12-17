#' @title POST function
#' @description Phylogeny-Guided Microbiome OTU-Specific Association Test
#' @param YY a numerical vector of outcome data with length n (sample size); can be binary outcome or continuous data
#' @param XX a matrix of covariant data with dimention n by q (number of covariants). Each column is a covariate.
#' @param OTU a matrix of OTU counts table with dimention n by p (number of OTUs).
#' @param tree a phylogenetic tree
#' @param cValues a vector of c values
#' @param trait character indicating the type of outcome, need to be "binomial" or "continuous".
#' @param verbose logical; generate progress screen prints
#' @author Caizhi Huang; Jung-Ying Tzeng
#' @import stats
#' @import MASS
#' @import CompQuadForm
#' @import ACAT
#' @references Debashis Ghosh and Xiang Zhan. "CKAT software" (2016)
#' @return Return a matrix with the rows corresponding to OTUs
#'   Columns correspond to: \cr
#'   (1) Name:       OTU name \cr
#'   (2) POST_rawp:  POST raw p-value; \cr
#'   (3) SO_rawp:    Single OTU test raw p-value\cr
#'   (4) minP:       The p-value for minC
#'   (5) minC:       The c value with minimal p-value among the c vectors \cr
#' @examples
#' if(!requireNamespace("GUniFrac",quietly = TRUE)){install.package("GUniFrac")}
#' if(!requireNamespace("GUniFrac",quietly = TRUE)){devtools::install_github("yaowuliu/ACAT")}
#' library("GUniFrac")
#' library("ACAT")
#' library("CompQuadForm")
#' data("throat.tree")
#' data("throat.otu.tab")
#' data("throat.meta")
#' YY <- (throat.meta[,"RespiratoryDiseaseStatus_severity_timeframe"]!="Healthy")+0
#' XX <- (throat.meta[,"SmokingStatus"]!="Smoker")+0
#' result <- POST(YY = YY, XX = XX, OTU = throat.otu.tab[,1:100], tree = throat.tree,
#'               cValues = seq(0,0.05,by=0.01), trait = 'binomial', verbose = TRUE)
#' @export


POST <- function( YY, #response
                  XX, #covariant
                  OTU, #OTU table
                  tree, #phylogeney three
                  cValues = seq(0,0.05,0.01),
                  trait = 'binomial', #continous or binomial
                  verbose = FALSE) {

  #options(warn = -1)
  options(scipen=999)
  
  OTU     <- as.matrix(OTU)
  nsubj   <- nrow(OTU)
  nOTU    <- ncol(OTU)

  XX <- as.matrix(XX)

  kernellist <- getkernellist(cValues = cValues,
                              tree    = tree,
                              OTU     = OTU,
                              verbose = verbose)


  postp <- numeric(length = nOTU)
  minp  <- numeric(length = nOTU)
  minc  <- numeric(length = nOTU)
  c0p   <- numeric(length = nOTU)

  for( m in 1L:nOTU) {

    if (verbose) cat("OTU", m, "\n")

    pvCKATv <- c()
    for (cv in 1L:length(cValues)) {

      if(trait == "binomial"){
        pvCKAT <- CKAT.bin(y=YY, K=kernellist[[cv]][[m]][[1]], X=XX)
        if(pvCKAT>1){pvCKAT=1}
      }
      else{
        pvCKAT <- CKAT.cont(y=YY, K=kernellist[[cv]][[m]][[1]], X=XX)
        if(pvCKAT>1){pvCKAT=1}
      }
      pvCKATv[cv] <- pvCKAT
    }

    if(sum(pvCKATv==1)*sum(pvCKATv==0)!=0){
      postp[m]=0
    }else{
        postp[m] <- ACAT(pvCKATv)
        }
    minp[m] <- min(pvCKATv)
    minc[m] <- cValues[which.min(pvCKATv)]
    c0p[m]  <- pvCKATv[1]
  }

  otuname                <- colnames(OTU)
  resultsMat             <- cbind(otuname, postp, c0p, minp, minc)
  Pvalue_Table           <- resultsMat[order(resultsMat[,"postp"]),]
  colnames(Pvalue_Table) <- c("Name", "POST_rawP", "SO_rawP", "minP", "minC")
  Pvalue_Table[,2:4]     <- matrix(format(as.numeric(Pvalue_Table[,2:4]), digits = 3,nsamll = 4,scientific = TRUE),ncol = 3)
  return(Pvalue_Table)
}
