#' mRSD checker
#'
#'
#'  
#' @param DataMatrix A matrix with feature groups as columns and samples for rows. 
#' @param subclasses A vector with characters or factors corresponding to those found in the 'labels' attribute. 
#' This can also be a list of vectors in case multiple variables are supplied in the labels parameter.
#' @param labels A vector of labels which can be used to group the samples together in subclasses. 
#' This attribute is necessary if 'subclasses' is provided. Note this can also be a list in case 2 variables need to 
#' be combined for the calculation (for example Batch x sample group). The length of 'labels', or each element in the list
#' has to be the same as the number of rows in the DataMatrix.
#' @param RSD_threshold The maximal RSD value (sd/mean) allowed in a single subclass.
#' @param ignore_missing_values Whether to ignore the missing values. If set to FALSE 
#' the calculation of the RSD will result in a lot of 'NAs' or 'NANs'. Only the features 
#' with no missing values will produce a numeric RSD.
#' @param NA_or_numeric_limit Define the missing values, this can be set to 'NA' or to a 
#' numerical value in which case all values below this value will be deemed missing.
#' @param feature_orientation Indicates whether the features can be found in the columns (default) or 
#' in the rows. With the default setting every row corresponds to a sample, and every column to a feature/group
#' @param groups_ok_threshold The amount of subclasses (as defined in subclasses_list)
#' @param include_RSD_matrix If set to 'TRUE' the output will not only be a vector with the indices
#' of groups that passed the RSD check, but also the matrix with all RSD values is provided. 
#'   
#' @return 
#' RSD_results A vector with the indices of the groups that passed the test for RSD. 
#' If however the 'include_RSD_matrix' attribute is set to 'TRUE' the results will be 
#' a list with the vector of indices (list entry 1) and the matrix of RSD for each 
#' combination of subclass and feature (list entry 2)       
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' 
#' @importFrom stats na.omit
#' 
#' @export
mRSD_check <- function(DataMatrix, subclasses = NULL,  labels = NULL, RSD_threshold = 0.4, ignore_missing_values = TRUE, NA_or_numeric_limit = NA, feature_orientation = "columns", groups_ok_threshold = 1, include_RSD_matrix = FALSE)    {
    
    
    # checks
    if(!feature_orientation %in% c("columns", "rows")){
        warning("'feature_orientation' was not one of two allowed possibilities. Set to the default 'columns'.")
        feature_orientation <- "columns"
    }
    
    if(feature_orientation == "rows"){
        DataMatrix <- t(DataMatrix)
    }
    
    MultiVar <- FALSE
    if(!is.null(labels)){
        if("list" %in% class(labels)){
            MultiVar <- TRUE
            labels.unique.list <- list()
            for(l in 1:length(labels)){
                if( length(labels[[l]]) != nrow(DataMatrix) ){
                    stop("the length of the labels vector is not equal to the amount of rows in the DataMatrix.")
                }
                labels.unique.list[[l]] <- unique(labels[[l]])
            }
        } else{
            if( length(labels) != nrow(DataMatrix) ){
                stop("the length of the labels vector is not equal to the amount of rows in the DataMatrix.")
            }
        }
    }
    
    if(is.null(labels) & is.null(subclasses)){
        # if the 'labels' and 'subclasses are both not supplied it is assumed that all rows belong to the same class
        subclasses_list <- list()
        subclasses_list[[1]] <- seq(1,nrow(DataMatrix))
    } else if(!is.null(labels) & is.null(subclasses)){
        subclasses_list <- list()
        if(MultiVar){
            combos <- expand.grid(labels.unique.list)
            for(k in 1:nrow(combos)){
                subset_true <- data.frame(matrix(FALSE, nrow = nrow(DataMatrix), ncol = ncol(combos)))
                for(l in 1:ncol(combos)){
                    subset_true[labels[[l]] == combos[k,l], l] = TRUE
                }
                subclasses_list[[k]] <- which(rowSums(subset_true) == ncol(combos))
            }
        } else{
            for(classes in 1:length(unique(labels))){
                subclasses_list[[classes]] <- which( labels == unique(labels)[classes])
            }
        }
    }
    
    if(!is.null(subclasses)){
        
        if(MultiVar){
            subclasses_list <- list()
            combos <- expand.grid(subclasses)
            for(k in 1:nrow(combos)){
                subset_true <- data.frame(matrix(FALSE, nrow = nrow(DataMatrix), ncol = ncol(combos)))
                for(l in 1:ncol(combos)){
                    subset_true[labels[[l]] == combos[k,l], l] = TRUE
                }
                subclasses_list[[k]] <- which(rowSums(subset_true) == ncol(combos))
            }
        } else{
            if( class(subclasses) %in% c("factor", "character") ){
                if(is.null(labels)){
                    stop("When submitting 'subclasses_list' with factors or characters the labels attribute should also be provided")
                }
                subclasses_list <- list()
                for(class_entries in 1:length(subclasses)){
                    subclasses_list[[class_entries]] <- which(labels %in% subclasses[class_entries])
                }
            } else{
                stop("Unexpected form for 'subclasses_list'. It is not a list of factors, characters, numerics or integers.")
            }
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
    N_subclasses_to_check <- length(subclasses_list)
    RSD_ok <-rep(0, ncol(DataMatrix))
    RSD_matrix <- matrix(0, ncol(DataMatrix), nrow = N_subclasses_to_check)
    for(subCL in 1:N_subclasses_to_check){
        
        subclass_Data <- DataMatrix[subclasses_list[[subCL]], ]
        
        if(nrow(subclass_Data) == 0){
            RSD_matrix[subCL,] <- 0
        } else{
            for(grp in 1:ncol(DataMatrix)){
                
                subclass_feat_Data <- as.numeric(subclass_Data[, grp])
                
                if(ignore_missing_values){
                    subclass_feat_Data <- na.omit(subclass_feat_Data) 
                }
                
                RSD_feat = sd(subclass_feat_Data) / mean(subclass_feat_Data)
                if(is.nan(RSD_feat) | is.na(RSD_feat) ){ RSD_feat <- -1}
                
                RSD_matrix[subCL,grp] <- RSD_feat
                if(RSD_feat >= 0 & RSD_feat <= RSD_threshold){
                    RSD_ok[grp] <- RSD_ok[grp] + 1
                } 
            }
        }
    }
    
    if(!include_RSD_matrix){
        RSD_results <- which(RSD_ok >= groups_ok_threshold)    
    } else{
        RSD_results <- list()
        RSD_results[[1]] <- which(RSD_ok >= groups_ok_threshold)
        RSD_results[[2]] <- RSD_matrix
        if(MultiVar){
            RSD_results[[3]] <- combos
        }
    }
    
    return(RSD_results)
    
    }
    