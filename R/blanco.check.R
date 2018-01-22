#' Blanco filter
#'
#' function only lets groups of samples/subclasses pass for which at least one of the 
#' subclasses has a median value x times (10 default) larger than the maximal blanco value
#'
#'  
#' @param DataMatrix A matrix with feature groups as columns and samples for rows. 
#' @param subclasses.list A list with either A) in each list entry a number of integers or 
#' numerics corresponding to the rows that form the subclasses to be checked, or B) in each 
#' list entry characters or factors corresponding to those found in the 'labels' attribute.
#' @param blanco.entry The label of the blanco samples (as to be found in the 'labels' attribute) or 
#' a vector of row numbers corresponding to the rows with blanco samples.
#' @param blanco.data If the blanco data is not in 'DataMatrix' but in an additional matrix than it can
#' be provided with this entry.
#' @param Blanco.multiplier The maximal RSD value (sd/mean) allowed in a single subclass.
#' @param ignore.missing.values Whether to ignore the missing values. If set to FALSE 
#' the calculation of the median intensities will result in a lot 'NAs' or 'NaNs'
#' @param labels A vector of labels which can be used to group the samples together in subclasses. 
#' This attribute is necessary if 'subclasses.list' consists of character or factor entries.
#' @param feature.orientation Indicates whether the features can be found in the columns (default) or 
#' in the rows. With the default setting every row corresponds to a sample, and every column to a feature/group
#' @param groups.ok.threshold The minimal amount of subclasses (as defined in subclasses.list) that have to pass the test
#' 
#'   
#' @return 
#' blanco.ok A vector with the indices of groups/features which passed the blanco filter.       
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' 
#' 
#' @export
#' 
#' @importFrom matrixStats colMaxs
#' @importFrom stats median na.omit
#' 
#' 
blanco.check = function(DataMatrix, subclasses.list = NULL, blanco.entry = NULL, blanco.data = NULL, Blanco.multiplier = 10,ignore.missing.values = TRUE, labels = NULL, feature.orientation = "columns", groups.ok.threshold = 1){
    
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
    
    if(is.null(blanco.entry) & is.null(blanco.data)){
        stop("Neither the 'blanco.entry' parameter nor the blanco.data entry is provided.")
    }
    
    if(!is.null(blanco.data)){
        if(ncol(blanco.data) != ncol(DataMatrix)){
            stop("The number of columns in 'blanco.data' is not equal to that of 'DataMatrix'.")
        }
        max.blancos <- matrixStats::colMaxs(blanco.data,na.rm=TRUE)
    } else if(is.character(blanco.entry) | is.factor(blanco.entry)){
        if(is.null(labels)){
            stop("The 'blanco.entry' is a character or factor but the 'labels' attribute is not provided.")
        }
        blanco.entry = which(labels %in% blanco.entry)
        max.blancos <- matrixStats::colMaxs(DataMatrix[blanco.entry,],na.rm=TRUE)
    }
    
    # Main function
    blanco.ok <-rep(0, ncol(DataMatrix))
    
    for(subCL in 1:N.sublasses.to.check){
        subclass.Data <- DataMatrix[subclasses.list[[subCL]], ]
        for(grp in 1:ncol(DataMatrix)){
            subclass.median <- median(na.omit(as.numeric(subclass.Data[, grp])))
            if(subclass.median >= Blanco.multiplier*max.blancos[grp]){
                blanco.ok[grp] <- blanco.ok[grp] + 1
            }
        }
    }
    
    blanco.ok <- which(blanco.ok >= groups.ok.threshold)
    return(blanco.ok)
    
}