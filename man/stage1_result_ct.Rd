\name{stage1_result_ct}
\alias{stage1_result_ct}
\docType{data}
\title{
The output from the function RunDeconvMethyl().
}
\description{
The function RunDeconvMethyl() returns a list of seven objects. The object stage1_result_ct is this file. These objects are the deconvoluted methylation profile of the samples, the proportions of each cell type, the iterations, the explained variance, the residual sum of squares, the Akaike information criterion (AIC), and the residual sum of squares per iteration in the deconvolution algorithm.}
\usage{data("stage1_result_ct")}
\format{
  The format is:
List of 7
 $ methylation       : num [1:499, 1:3] 0.1607 0.8006 0.6967 0.524 0.0751 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:499] "cg00323915" "cg21830221" "cg11348106" "cg12232731" ...
  .. ..$ : NULL
 $ proportions       : num [1:80, 1:3] 0.31 0.55 0.86 0.58 0.73 0.14 0.49 0 0 0.67 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:80] "TCGA-OU-A5PI-01A" "TCGA-OR-A5K8-01A" "TCGA-OR-A5JX-01A" "TCGA-OR-A5LT-01A" ...
  .. ..$ : NULL
 $ iterations        : num 234
 $ explained.variance: num 0.955
 $ res.sum.squares   : num 1039
 $ aic               : num -264674
 $ rss.per.iteration : num [1:234] 1366 1165 1098 1071 1061 ...
}
\examples{
data(stage1_result_ct)
## maybe str(stage1_result_ct) ; plot(stage1_result_ct) ...
}
\keyword{datasets}
