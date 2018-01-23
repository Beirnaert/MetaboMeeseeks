#' Missing value proportion checker
#'
#' function only lets groups of samples/subclasses pass which have a proportion of
#' missing values less than the allowed threshold
#'
#' @param DataMatrix A matrix with feature groups as columns and samples for rows. 
#' @param subclasses.list A list with either A) in each list entry a number of integers or 
#' numerics corresponding to the rows that form the subclasses to be checked, or B) in each 
#' list entry characters or factors corresponding to those found in the 'labels' attribute.
#' @param prop.missing.threshold The maximal proportion of missing values allowed in a single subclass.
#' @param Na.or.numeric.limit Define the missing values, this can be set to 'NA' or to a 
#' numerical value in which case all values below this value will be deemed missing.
#' @param labels A vector of labels which can be used to group the samples together in subclasses. 
#' This attribute is necessary if 'subclasses.list' consists of character or factor entries.
#' @param feature.orientation Indicates whether the features can be found in the 
#' columns (default) or in the rows. With the default setting every row corresponds 
#' to a sample, and every column to a feature/group
#' @param groups.ok.threshold The amount of subclasses (as defined in subclasses.list)
#' 
#' @return 
#' Groups.that.passed A vector with the indices of the groups that passed the test 
#' for proportion of missing values.
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#'
#' 
#' @export
missing.values.check = function(DataMatrix, subclasses.list = NULL, prop.missing.threshold = 0.2, Na.or.numeric.limit = NA, labels = NULL, feature.orientation = "columns", groups.ok.threshold = 1){

   
    # checks
    if(!feature.orientation %in% c("columns", "rows")){
        warning("'feature.orientation' was not one of two allowed possibilities. Set to the default 'columns'.")
        feature.orientation = "columns"
    }
    
    if(feature.orientation == "rows"){
        DataMatrix = t(DataMatrix)
    }
    
    if(!is.null(labels)){
        if( length(labels) != nrow(DataMatrix) ){
            stop("the length of the labels vector is not equal to the amount of rows in the DataMatrix.")
        }
    }
    
    if(is.null(labels) & is.null(subclasses.list)){
        # if the 'labels' and 'subclasses.list' are both not supplied it is assumed that all raws belong to the same class
        subclasses.list = list()
        subclasses.list[[1]] <- seq(1,nrow(DataMatrix))
    } else if(!is.null(labels) & is.null(subclasses.list)){
        subclasses.list <- list()
        for(classes in 1:length(unique(labels))){
            subclasses.list[[classes]] <- unique(labels)[classes]
        }
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
    
    if(!is.na(Na.or.numeric.limit) & !is.numeric(Na.or.numeric.limit)){
        warning("The 'Na.or.numeric.limit' attribute is neither 'NA' nor a numeric value. Set to the defualt 'NA'.")
        Na.or.numeric.limit <- NA
    } 
    
    if(is.numeric(Na.or.numeric.limit)){
        DataMatrix[DataMatrix<Na.or.numeric.limit] <- NA
    }
    
   
    # Main function
    
    subgroups.ok <- rep(0,ncol(DataMatrix))
    
    for(subCL in 1:N.sublasses.to.check){
        
        subclass.Data <- DataMatrix[subclasses.list[[subCL]], ]
        
        for(grp in 1:ncol(DataMatrix)){
            if(sum(is.na(subclass.Data[, grp]))/nrow(subclass.Data) <= prop.missing.threshold ){
                subgroups.ok[grp] <- subgroups.ok[grp] + 1
            }
        }
        
    }
    
    
    Groups.that.passed = which(subgroups.ok >= groups.ok.threshold)
    
    return(Groups.that.passed)
}