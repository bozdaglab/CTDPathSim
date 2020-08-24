\name{FindSim}
\alias{FindSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to find similarity score between patient and cell line.
}
\description{
FindSim uses the methylation profile, expression profile, and pathway for both the patient and cell line to compute a similarity score.
}
\usage{
FindSim(pat_de_genes,pat_reactome,expression,cell_de_genes,cell_reactome,ccle_expr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{pat_de_genes}{
A R dataframe object with all of the patients DE genes.
    }\item{pat_reactome}{
A R dataframe object with the patients reactome pathways.
    }\item{expression}{
A R dataframe object with the patients expression values for each gene and cell type.
    }\item{cell_de_genes}{
A R dataframe object with all of the cel line's DE genes.
    }\item{cell_reactome}{
A R dataframe object with the cell line's reactome pathways.
    }\item{ccle_expr}{
A R dataframe object with the cell line's expression values for each gene and cell type.
    }
}
\value{
A R Dataframe object that includes the similarity score for the patient and cell line
}
\author{
Banabithi Bose
}
\examples{
load("~/CTDPathSim/inst/extdata/Example.patient1.de.genes.rda")
load("~/CTDPathSim/inst/extdata/Example.patient1.pathways.rda")
load("~/CTDPathSim/inst/extdata/Example.patient1.deconvoluted.expression.rda")
load("~/CTDPathSim/inst/extdata/Example.cellLine1.de.genes.rda")
load("~/CTDPathSim/inst/extdata/Example.cellLine1.pathways.rda")
load("~/CTDPathSim/inst/extdata/Example.cellLine1.expression.rda")
patient1.cellLine1.spearman<-FindSim(pat_de_genes,pat_reactome,expression,
cell_de_genes,cell_reactome,ccle_expr)
}
