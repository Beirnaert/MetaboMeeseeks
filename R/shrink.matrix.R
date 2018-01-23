#' Shrink matrix from multiple samples per class to 1 value per class
#'
#' This function takes a matrix with with multiple samples per subclass and converts these to 
#' mean (or other metric) and sd values per class. Effectively it is just a wrapper for a few functions
#' of the matrixStats package to combine the results for multiple classes. I.e. the same can easily be
#' achieved by for instance rbinding two results from colMedians with two different row subsets
#'
#'  
#' @param DataMatrix A matrix with feature groups as columns and samples for rows. 
#' @param subclasses.list A list with either A) in each list entry a number of integers or 
#' numerics corresponding to the rows that form the subclasses to be checked, or B) in each 
#' list entry characters or factors corresponding to those found in the 'labels' attribute.
#' @param ignore.missing.values Whether to ignore the missing values. If set to FALSE 
#' the calculation of the median intensities will result in a lot 'NAs' or 'NaNs'
#' @param labels A vector of labels which can be used to group the samples together in subclasses. 
#' This attribute is necessary if 'subclasses.list' consists of character or factor entries.
#' @param feature.orientation Indicates whether the features can be found in the columns (default) or 
#' in the rows. With the default setting every row corresponds to a sample, and every column to a feature/group
#' @param metric The metric to be used, can be any of 'mean', 'median', 'max', 'min', 'var' or 'sd'.
#' 
#'   
#' @return 
#' Shrunken.matrix The shrunken matrix       
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' #Generate random matrix
#' testmatrix <- matrix(runif(n=90), nrow = 9, ncol = 10)
#' sample.groups <- c(rep("A",3), rep("B",3), rep("C",3))
#' subclasses <- as.list(unique(c(rep("A",3), rep("B",3), rep("C",3))))
#' 
#' Shrunken.matrix = shrink.matrix(DataMatrix = testmatrix, 
#'                                 subclasses.list = subclasses, 
#'                                 labels = sample.groups, 
#'                                 metric = "median")
#' 
#' #Note that the exact same result can be obtained with the following
#' library(matrixStats)
#' Shrunken.matrix2 = rbind(matrixStats::colMedians(testmatrix, 
#'                                                  rows = which(sample.groups == subclasses[[1]]), 
#'                                                  na.rm = TRUE), 
#'                          matrixStats::colMedians(testmatrix, 
#'                                                  rows = which(sample.groups == subclasses[[2]]), 
#'                                                  na.rm = TRUE),
#'                          matrixStats::colMedians(testmatrix, 
#'                                                  rows = which(sample.groups == subclasses[[3]]), 
#'                                                  na.rm = TRUE))
#' 
#' @export
#' 
#' @import matrixStats
#' 
shrink.matrix = function(DataMatrix, subclasses.list = NULL, ignore.missing.values = TRUE, labels = NULL, feature.orientation = "columns", metric = "median"){
    
    # checks
    if(!feature.orientation %in% c("columns", "rows")){
        warning("'feature.orientation' was not one of two allowed possibilities. Set to the default 'columns'.")
        feature.orientation = "columns"
    }
    
    if(feature.orientation == "rows"){
        DataMatrix = t(DataMatrix)
        if(!is.null(blanco.data)){
            blanco.data = t(blanco.data)
        }
    }
    
    if(!is.null(labels)){
        if( length(labels) != nrow(DataMatrix) ){
            stop("the length of the labels vector is not equal to the amount of rows in the DataMatrix.")
        }
    }
    
    if(is.null(labels) & is.null(subclasses.list)){
        stop("Neither 'labels' nor 'subclasses.list' is provided.")
    }
    
    if(!is.null(subclasses.list)){
        if( class(subclasses.list[[1]][1]) %in% c("factor", "character") ){
            if(is.null(labels)){
                stop("When submitting 'subclasses.list' with factors or characters the labels attribute should also be provided")
            }
            subclasses.list.numeric = list()
            for(list.entries in 1:length(subclasses.list)){
                subclasses.list.numeric[[list.entries]] <- which(labels %in% subclasses.list[[list.entries]])
            }
            subclasses.list <- subclasses.list.numeric
            N.sublasses.to.check <- length(subclasses.list)
        } else if( class(subclasses.list[[1]][1]) %in% c("numeric", "integer")){
            N.sublasses.to.check <- length(subclasses.list)
        } else{
            stop("Unexpected form for 'subclasses.list'. It is not a list of factors, characters, numerics or integers.")
        }
    }
    
    
    # Main function
    Shrunken.matrix <-matrix(0, nrow = N.sublasses.to.check, ncol = ncol(DataMatrix))
    
    for(subCL in 1:N.sublasses.to.check){
        
        if(metric == "median"){
            Shrunken.matrix[subCL,] = matrixStats::colMedians(DataMatrix, rows = subclasses.list[[subCL]], na.rm = ignore.missing.values)
        }
        if(metric == "mean"){
            Shrunken.matrix[subCL,] = colMeans(DataMatrix[subclasses.list[[subCL]],], na.rm = ignore.missing.values)
        }
        if(metric == "min"){
            Shrunken.matrix[subCL,] = matrixStats::colMins(DataMatrix, rows = subclasses.list[[subCL]], na.rm = ignore.missing.values)
        }
        if(metric == "max"){
            Shrunken.matrix[subCL,] = matrixStats::colMaxs(DataMatrix, rows = subclasses.list[[subCL]], na.rm = ignore.missing.values)
        }
        if(metric == "var" | metric == "variance"){
            Shrunken.matrix[subCL,] = matrixStats::colVars(DataMatrix, rows = subclasses.list[[subCL]], na.rm = ignore.missing.values)
        }
        if(metric == "sd" | metric == "standard deviation"){
            Shrunken.matrix[subCL,] = matrixStats::colSds(DataMatrix, rows = subclasses.list[[subCL]], na.rm = ignore.missing.values)
        }
    }
    
    return(Shrunken.matrix)
    
}
