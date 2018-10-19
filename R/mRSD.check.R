#' mRSD checker
#'
#' function only lets groups of samples/subclasses pass which have a proportion of
#' missing values less than the allowed threshold
#'
#'  
#' @param DataMatrix A matrix with feature groups as columns and samples for rows. 
#' @param subclasses A vector with characters or factors corresponding to those found in the 'labels' attribute. 
#' This can also be a list of vectors in case multiple variables are supplied in the labels parameter.
#' @param labels A vector of labels which can be used to group the samples together in subclasses. 
#' This attribute is necessary if 'subclasses' is provided. Note this can also be a list in case 2 variables need to 
#' be combined for the calculation (for example Batch x sample group). The length of 'labels', or each element in the list
#' has to be the same as the number of rows in the DataMatrix.
#' @param RSD.threshold The maximal RSD value (sd/mean) allowed in a single subclass.
#' @param ignore.missing.values Whether to ignore the missing values. If set to FALSE 
#' the calculation of the RSD will result in a lot of 'NAs' or 'NANs'. Only the features 
#' with no missing values will produce a numeric RSD.
#' @param Na.or.numeric.limit Define the missing values, this can be set to 'NA' or to a 
#' numerical value in which case all values below this value will be deemed missing.
#' @param feature.orientation Indicates whether the features can be found in the columns (default) or 
#' in the rows. With the default setting every row corresponds to a sample, and every column to a feature/group
#' @param groups.ok.threshold The amount of subclasses (as defined in subclasses.list)
#' @param include.RSD.matrix If set to 'TRUE' the output will not only be a vector with the indices
#' of groups that passed the RSD check, but also the matrix with all RSD values is provided. 
#'   
#' @return 
#' RSD.results A vector with the indices of the groups that passed the test for RSD. 
#' If however the 'include.RSD.matrix' attribute is set to 'TRUE' the results will be 
#' a list with the vector of indices (list entry 1) and the matrix of RSD for each 
#' combination of subclass and feature (list entry 2)       
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' 
#' @importFrom stats na.omit
#' 
#' @export
mRSD.check = function(DataMatrix, subclasses = NULL,  labels = NULL, RSD.threshold = 0.4,ignore.missing.values = TRUE, Na.or.numeric.limit = NA, feature.orientation = "columns", groups.ok.threshold = 1, include.RSD.matrix = FALSE)    {
    
    
    # checks
    if(!feature.orientation %in% c("columns", "rows")){
        warning("'feature.orientation' was not one of two allowed possibilities. Set to the default 'columns'.")
        feature.orientation = "columns"
    }
    
    if(feature.orientation == "rows"){
        DataMatrix = t(DataMatrix)
    }
    
    MultiVar = FALSE
    if(!is.null(labels)){
        if("list" %in% class(labels)){
            MultiVar = TRUE
            labels.unique.list = list()
            for(l in 1:length(labels)){
                if( length(labels[[l]]) != nrow(DataMatrix) ){
                    stop("the length of the labels vector is not equal to the amount of rows in the DataMatrix.")
                }
                labels.unique.list[[l]] = unique(labels[[l]])
            }
        } else{
            if( length(labels) != nrow(DataMatrix) ){
                stop("the length of the labels vector is not equal to the amount of rows in the DataMatrix.")
            }
        }
    }
    
    if(is.null(labels) & is.null(subclasses)){
        # if the 'labels' and 'subclasses are both not supplied it is assumed that all rows belong to the same class
        subclasses.list = list()
        subclasses.list[[1]] <- seq(1,nrow(DataMatrix))
    } else if(!is.null(labels) & is.null(subclasses)){
        subclasses.list <- list()
        if(MultiVar){
            combos = expand.grid(labels.unique.list)
            for(k in 1:nrow(combos)){
                subset.true = data.frame(matrix(FALSE, nrow = nrow(DataMatrix), ncol = ncol(combos)))
                for(l in 1:ncol(combos)){
                    subset.true[labels[[l]] == combos[k,l], l] = TRUE
                }
                subclasses.list[[k]] = which(rowSums(subset.true) == ncol(combos))
            }
        } else{
            for(classes in 1:length(unique(labels))){
                subclasses.list[[classes]] <- which( labels == unique(labels)[classes])
            }
        }
    }
    
    if(!is.null(subclasses)){
        
        if(MultiVar){
            subclasses.list <- list()
            combos = expand.grid(subclasses)
            for(k in 1:nrow(combos)){
                subset.true = data.frame(matrix(FALSE, nrow = nrow(DataMatrix), ncol = ncol(combos)))
                for(l in 1:ncol(combos)){
                    subset.true[labels[[l]] == combos[k,l], l] = TRUE
                }
                subclasses.list[[k]] = which(rowSums(subset.true) == ncol(combos))
            }
        } else{
            if( class(subclasses) %in% c("factor", "character") ){
                if(is.null(labels)){
                    stop("When submitting 'subclasses.list' with factors or characters the labels attribute should also be provided")
                }
                subclasses.list = list()
                for(class.entries in 1:length(subclasses)){
                    subclasses.list[[class.entries]] <- which(labels %in% subclasses[class.entries])
                }
            } else{
                stop("Unexpected form for 'subclasses.list'. It is not a list of factors, characters, numerics or integers.")
            }
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
    N.subclasses.to.check <- length(subclasses.list)
    RSD.ok <-rep(0, ncol(DataMatrix))
    RSD.matrix <- matrix(0, ncol(DataMatrix), nrow = N.subclasses.to.check)
    for(subCL in 1:N.subclasses.to.check){
        
        subclass.Data <- DataMatrix[subclasses.list[[subCL]], ]
        
        if(nrow(subclass.Data) == 0){
            RSD.matrix[subCL,] = 0
        } else{
            for(grp in 1:ncol(DataMatrix)){
                
                subclass.feat.Data <- as.numeric(subclass.Data[, grp])
                
                if(ignore.missing.values){
                    subclass.feat.Data <- na.omit(subclass.feat.Data) 
                }
                
                RSD.feat = sd(subclass.feat.Data) / mean(subclass.feat.Data)
                if(is.nan(RSD.feat) | is.na(RSD.feat) ){ RSD.feat <- -1}
                
                RSD.matrix[subCL,grp] <- RSD.feat
                if(RSD.feat >= 0 & RSD.feat <= RSD.threshold){
                    RSD.ok[grp] <- RSD.ok[grp] + 1
                } 
            }
        }
    }
    
    if(!include.RSD.matrix){
        RSD.results <- which(RSD.ok >= groups.ok.threshold)    
    } else{
        RSD.results <- list()
        RSD.results[[1]] <- which(RSD.ok >= groups.ok.threshold)
        RSD.results[[2]] <- RSD.matrix
        if(MultiVar){
            RSD.results[[3]] <- combos
        }
    }
    
    return(RSD.results)
    
    }
    