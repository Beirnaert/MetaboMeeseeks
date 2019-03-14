#' Missing value proportion checker
#'
#' function only lets groups of samples/subclasses pass which have a proportion of
#' missing values less than the allowed threshold
#'
#' @param DataMatrix A matrix with feature groups as columns and samples for rows. 
#' @param subclasses_list A list with either A) in each list entry a number of integers or 
#' numerics corresponding to the rows that form the subclasses to be checked, or B) in each 
#' list entry characters or factors corresponding to those found in the 'labels' attribute.
#' @param prop_missing_threshold The maximal proportion of missing values allowed in a single subclass.
#' @param NA_or_numeric_limit Define the missing values, this can be set to 'NA' or to a 
#' numerical value in which case all values below this value will be deemed missing.
#' @param labels A vector of labels which can be used to group the samples together in subclasses. 
#' This attribute is necessary if 'subclasses_list' consists of character or factor entries.
#' @param feature_orientation Indicates whether the features can be found in the 
#' columns (default) or in the rows. With the default setting every row corresponds 
#' to a sample, and every column to a feature/group
#' @param groups_ok_threshold The amount of subclasses (as defined in subclasses_list)
#' 
#' @return 
#' Groups_that_passed A vector with the indices of the groups that passed the test 
#' for proportion of missing values.
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#'
#' 
#' @export
MissingValuesCheckR = function(DataMatrix, subclasses_list = NULL, prop_missing_threshold = 0.2, NA_or_numeric_limit = NA, labels = NULL, feature_orientation = "columns", groups_ok_threshold = 1){

   
    # checks
    if(!feature_orientation %in% c("columns", "rows")){
        warning("'feature_orientation' was not one of two allowed possibilities. Set to the default 'columns'.")
        feature_orientation = "columns"
    }
    
    if(feature_orientation == "rows"){
        DataMatrix = t(DataMatrix)
    }
    
    if(!is.null(labels)){
        if( length(labels) != nrow(DataMatrix) ){
            stop("the length of the labels vector is not equal to the amount of rows in the DataMatrix.")
        }
    }
    
    if(is.null(labels) & is.null(subclasses_list)){
        # if the 'labels' and 'subclasses_list' are both not supplied it is assumed that all raws belong to the same class
        subclasses_list = list()
        subclasses_list[[1]] <- seq(1,nrow(DataMatrix))
    } else if(!is.null(labels) & is.null(subclasses_list)){
        subclasses_list <- list()
        for(classes in 1:length(unique(labels))){
            subclasses_list[[classes]] <- unique(labels)[classes]
        }
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
            N.sublasses.to.check <- length(subclasses_list)
        } else if( class(subclasses_list[[1]][1]) %in% c("numeric", "integer")){
            N.sublasses.to.check <- length(subclasses_list)
        } else{
            stop("Unexpected form for 'subclasses_list'. It is not a list of factors, characters, numerics or integers.")
        }
    }
    
    if(!is.na(NA_or_numeric_limit) & !is.numeric(NA_or_numeric_limit)){
        warning("The 'NA_or_numeric_limit' attribute is neither 'NA' nor a numeric value. Set to the defualt 'NA'.")
        NA_or_numeric_limit <- NA
    } 
    
    if(is.numeric(NA_or_numeric_limit)){
        DataMatrix[DataMatrix<NA_or_numeric_limit] <- NA
    }
    
   
    # Main function
    
    subgroups_ok <- rep(0,ncol(DataMatrix))
    
    for(subCL in 1:N.sublasses.to.check){
        
        subclass_Data <- DataMatrix[subclasses_list[[subCL]], ]
        
        for(grp in 1:ncol(DataMatrix)){
            if(sum(is.na(subclass_Data[, grp]))/nrow(subclass_Data) <= prop_missing_threshold ){
                subgroups_ok[grp] <- subgroups_ok[grp] + 1
            }
        }
        
    }
    
    
    Groups_that_passed = which(subgroups_ok >= groups_ok_threshold)
    
    return(Groups_that_passed)
}