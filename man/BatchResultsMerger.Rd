% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BatchResultsMerger.R
\name{BatchResultsMerger}
\alias{BatchResultsMerger}
\title{Restructure BatchProcessing and calculate combined results}
\usage{
BatchResultsMerger(PLSDAresults.list, BatchProcessingStructure,
  original.column.indices.list, GroupData.list, FoldChanges.list)
}
\arguments{
\item{PLSDAresults.list}{list of PLSDA variable importances}

\item{BatchProcessingStructure}{A BatchProcessingStructure}

\item{original.column.indices.list}{list with the original column indices}

\item{GroupData.list}{list with the original GroupData}

\item{FoldChanges.list}{list with fold change vectors}
}
\value{
BatchProcessingStructure A data frame with the list structure (quantifies the list nesting)
}
\description{
Restructure BatchProcessing and calculate combined results
}
\author{
Charlie Beirnaert, \email{charlie.beirnaert@uantwerpen.be}
}
