#' Blanco filter
#'
#' function only lets groups of samples/subclasses pass for which at least one of the 
#' subclasses has a median value x times (10 default) larger than the maximal blanco value
#'
#'  
#' @param DataMatrix A matrix with feature groups as columns and samples for rows. 
#' @param subclasses_list A list with either A) in each list entry a number of integers or 
#' numerics corresponding to the rows that form the subclasses to be checked, or B) in each 
#' list entry characters or factors corresponding to those found in the 'labels' attribute.
#' @param blank_entry The label of the blanco samples (as to be found in the 'labels' attribute) or 
#' a vector of row numbers corresponding to the rows with blanco samples.
#' @param blank_data If the blanco data is not in 'DataMatrix' but in an additional matrix than it can
#' be provided with this entry.
#' @param blank_multiplier The maximal RSD value (sd/mean) allowed in a single subclass.
#' @param ignore_missing_values Whether to ignore the missing values. If set to FALSE 
#' the calculation of the median intensities will result in a lot 'NAs' or 'NaNs'
#' @param labels A vector of labels which can be used to group the samples together in subclasses. 
#' This attribute is necessary if 'subclasses_list' consists of character or factor entries.
#' @param feature_orientation Indicates whether the features can be found in the columns (default) or 
#' in the rows. With the default setting every row corresponds to a sample, and every column to a feature/group
#' @param groups_ok_threshold The minimal amount of subclasses (as defined in subclasses_list) that have to pass the test
#' 
#'   
#' @return 
#' blank_ok A vector with the indices of groups/features which passed the blanco filter.       
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' 
#' @export
#' 
#' @importFrom matrixStats colMaxs
#' @importFrom stats median na.omit
#' 
#' 
BlankCheckR <- function(DataMatrix, subclasses_list = NULL, blank_entry = NULL, blank_data = NULL, blank_multiplier = 10,ignore_missing_values = TRUE, labels = NULL, feature_orientation = "columns", groups_ok_threshold = 1){
    
    # checks
    if(!feature_orientation %in% c("columns", "rows")){
        warning("'feature_orientation' was not one of two allowed possibilities. Set to the default 'columns'.")
        feature_orientation <- "columns"
    }
    
    if(feature_orientation == "rows"){
        DataMatrix <- t(DataMatrix)
        if(!is.null(blank_data)){
            blank_data <- t(blank_data)
        }
    }
    
    if(!is.null(labels)){
        if( length(labels) != nrow(DataMatrix) ){
            stop("the length of the labels vector is not equal to the amount of rows in the DataMatrix.")
        }
    }
    
    if(is.null(labels) & is.null(subclasses_list)){
        stop("Neither 'labels' nor 'subclasses_list' is provided.")
    }
    
    if(!is.null(subclasses_list)){
        if( class(subclasses_list[[1]][1]) %in% c("factor", "character") ){
            if(is.null(labels)){
                stop("When submitting 'subclasses_list' with factors or characters the labels attribute should also be provided")
            }
            subclasses_list_numeric = list()
            for(list_entries in 1:length(subclasses_list)){
                subclasses_list_numeric[[list_entries]] <- which(labels %in% subclasses_list[[list_entries]])
            }
            subclasses_list <- subclasses_list_numeric
            N_sublasses_to_check <- length(subclasses_list)
        } else if( class(subclasses_list[[1]][1]) %in% c("numeric", "integer")){
            N_sublasses_to_check <- length(subclasses_list)
        } else{
            stop("Unexpected form for 'subclasses_list'. It is not a list of factors, characters, numerics or integers.")
        }
    }
    
    if(is.null(blank_entry) & is.null(blank_data)){
        stop("Neither the 'blank_entry' parameter nor the blank_data entry is provided.")
    }
    
    if(!is.null(blank_data)){
        if(ncol(blank_data) != ncol(DataMatrix)){
            stop("The number of columns in 'blank_data' is not equal to that of 'DataMatrix'.")
        }
        max_blanks <- matrixStats::colMaxs(blank_data, na.rm=TRUE)
    } else if(is.character(blank_entry) | is.factor(blank_entry)){
        if(is.null(labels)){
            stop("The 'blank_entry' is a character or factor but the 'labels' attribute is not provided.")
        }
        blank_entry <- which(labels %in% blank_entry)
        max_blanks <- matrixStats::colMaxs(DataMatrix[blank_entry,],na.rm=TRUE)
    }
    
    # Main function
    blank_ok <-rep(0, ncol(DataMatrix))
    
    for(subCL in 1:N_sublasses_to_check){
        subclass_Data <- DataMatrix[subclasses_list[[subCL]], ]
        for(grp in 1:ncol(DataMatrix)){
            subclass_median <- median(na.omit(as.numeric(subclass_Data[, grp])))
            if(subclass_median >= blank_multiplier*max_blanks[grp]){
                blank_ok[grp] <- blank_ok[grp] + 1
            }
        }
    }
    
    blank_ok <- which(blank_ok >= groups_ok_threshold)
    return(blank_ok)
    
}