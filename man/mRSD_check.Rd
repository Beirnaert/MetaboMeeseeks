% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mRSD_check.R
\name{mRSD_check}
\alias{mRSD_check}
\title{mRSD checker}
\usage{
mRSD_check(DataMatrix, subclasses = NULL, labels = NULL,
  RSD_threshold = 0.4, ignore_missing_values = TRUE,
  NA_or_numeric_limit = NA, feature_orientation = "columns",
  groups_ok_threshold = 1, include_RSD_matrix = FALSE)
}
\arguments{
\item{DataMatrix}{A matrix with feature groups as columns and samples for rows.}

\item{subclasses}{A vector with characters or factors corresponding to those found in the 'labels' attribute. 
This can also be a list of vectors in case multiple variables are supplied in the labels parameter.}

\item{labels}{A vector of labels which can be used to group the samples together in subclasses. 
This attribute is necessary if 'subclasses' is provided. Note this can also be a list in case 2 variables need to 
be combined for the calculation (for example Batch x sample group). The length of 'labels', or each element in the list
has to be the same as the number of rows in the DataMatrix.}

\item{RSD_threshold}{The maximal RSD value (sd/mean) allowed in a single subclass.}

\item{ignore_missing_values}{Whether to ignore the missing values. If set to FALSE 
the calculation of the RSD will result in a lot of 'NAs' or 'NANs'. Only the features 
with no missing values will produce a numeric RSD.}

\item{NA_or_numeric_limit}{Define the missing values, this can be set to 'NA' or to a 
numerical value in which case all values below this value will be deemed missing.}

\item{feature_orientation}{Indicates whether the features can be found in the columns (default) or 
in the rows. With the default setting every row corresponds to a sample, and every column to a feature/group}

\item{groups_ok_threshold}{The amount of subclasses (as defined in subclasses_list)}

\item{include_RSD_matrix}{If set to 'TRUE' the output will not only be a vector with the indices
of groups that passed the RSD check, but also the matrix with all RSD values is provided.}
}
\value{
RSD_results A vector with the indices of the groups that passed the test for RSD. 
If however the 'include_RSD_matrix' attribute is set to 'TRUE' the results will be 
a list with the vector of indices (list entry 1) and the matrix of RSD for each 
combination of subclass and feature (list entry 2)
}
\description{
mRSD checker
}
\author{
Charlie Beirnaert, \email{charlie.beirnaert@uantwerpen.be}
}
