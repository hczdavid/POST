% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POST.R
\name{POST}
\alias{POST}
\title{POST function}
\usage{
POST(
  YY,
  XX,
  OTU,
  tree,
  cValues = seq(0, 0.05, 0.01),
  trait = "binomial",
  verbose = FALSE
)
}
\arguments{
\item{YY}{a numerical vector of outcome data with length n (sample size); can be binary outcome or continuous data}

\item{XX}{a matrix of covariant data with dimention n by q (number of covariants). Each column is a covariate.}

\item{OTU}{a matrix of OTU counts table with dimention n by p (number of OTUs).}

\item{tree}{a phylogenetic tree}

\item{cValues}{a vector of c values}

\item{trait}{character indicating the type of outcome, need to be "binomial" or "continuous".}

\item{verbose}{logical; generate progress screen prints}
}
\value{
Return a matrix with the rows corresponding to OTUs
  Columns correspond to: \cr
  (1) Name:       OTU name \cr
  (2) POST_rawp:  POST raw p-value; \cr
  (3) SO_rawp:    Single OTU test raw p-value\cr
  (4) minP:       The p-value for minC
  (5) minC:       The c value with minimal p-value among the c vectors \cr
}
\description{
Phylogeny-Guided Microbiome OTU-Specific Association Test
}
\examples{
if(!requireNamespace("GUniFrac",quietly = TRUE)){install.package("GUniFrac")}
if(!requireNamespace("GUniFrac",quietly = TRUE)){devtools::install_github("yaowuliu/ACAT")}
library("GUniFrac")
library("ACAT")
library("CompQuadForm")
data("throat.tree")
data("throat.otu.tab")
data("throat.meta")
YY <- (throat.meta[,"RespiratoryDiseaseStatus_severity_timeframe"]!="Healthy")+0
XX <- (throat.meta[,"SmokingStatus"]!="Smoker")+0
result <- POST(YY = YY, XX = XX, OTU = throat.otu.tab[,1:100], tree = throat.tree,
              cValues = seq(0,0.05,by=0.01), trait = 'binomial', verbose = TRUE)
}
\references{
Debashis Ghosh and Xiang Zhan. "CKAT software" (2016)
}
\author{
Caizhi Huang; Jung-Ying Tzeng
}
