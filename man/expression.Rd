\name{expression}
\alias{expression}
\docType{data}
\title{
R dataframe object with expression profile}
\description{
Example Expression profile
}
\usage{data("expression")}
\format{
  The format is:
 num [1:37109, 1:3] 0.1046 0.00333 0.34627 0.02466 0.01724 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:37109] "TSPAN6" "TNMD" "DPM1" "SCYL3" ...
  ..$ : chr [1:3] "Monocytes" "Vascular_endothelial_cells" "Adipocytes"
}
\source{
TCGA Database (https://portal.gdc.cancer.gov/)
}
xamples{
data(expression)
## maybe str(expression) ; plot(expression) ...
}
\keyword{datasets}
