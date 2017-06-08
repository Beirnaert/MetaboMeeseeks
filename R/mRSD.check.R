#' mRSD checker
#'
#' function only lets groups of samples/subclasses pass which have a proportion of
#' missing values less than the allowed threshold
#'
#'  
#' @param DataMatrix A matrix with feature groups as columns and samples for rows. 
#' @param subclasses.list A list with either A) in each list entry a number of integers or 
#' numerics corresponding to the rows that form the subclasses to be checked, or B) in each 
#' list entry characters or factors corresponding to those found in the 'labels' attribute.
#' @param RSD.threshold The maximal RSD value (sd/mean) allowed in a single subclass.
#' @param ignore.missing.values Whether to ignore the missing values. If set to FALSE 
#' the calculation of the RSD will result in a lot of 'NAs' or 'NANs'. Only the features 
#' with no missing values will produce a numeric RSD.
#' @param Na.or.numeric.limit Define the missing values, this can be set to 'NA' or to a 
#' numerical value in which case all values below this value will be deemed missing.
#' @param labels A vector of labels which can be used to group the samples together in subclasses. 
#' This attribute is necessary if 'subclasses.list' consists of character or factor entries.
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
#' @examples
#' # this is an example for 3 meaningless spectra
#' lengths_of_spectra <- c(100,150,120)
#' measurement_distance <- 0.01
#' starting_ppm_values <- c(8.7, 9.0, 9.0)
#' spectra <- list()
#' ppm_values <- list()
#' for (k in 1:3) {
#'     spectra[[k]] <- runif(lengths_of_spectra[k], min = 0, max = 10)
#'     
#'     # note the minus sign in the 'by' statement
#'     ppm_values[[k]] <- seq(from = starting_ppm_values[k], by = -measurement_distance, 
#'                            length.out = lengths_of_spectra[k])  
#' }
#' new.Data <- BuildRawDataMatrix(spectrum.list = spectra, ppm.list = ppm_values)
#' spectraMatrix <- new.Data$DataMatrix
#' ppmMatrix <- new.Data$ppmMatrix
#' 
#' @export
mRSD.check = function(DataMatrix, subclasses.list = NULL, RSD.threshold = 0.4,ignore.missing.values = TRUE, Na.or.numeric.limit = NA, labels = NULL, feature.orientation = "columns", groups.ok.threshold = 1, include.RSD.matrix = FALSE)    {
    
    
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
    
    RSD.ok <-rep(0, ncol(DataMatrix))
    RSD.matrix <- matrix(0, ncol(DataMatrix), nrow = N.sublasses.to.check)
    
    for(subCL in 1:N.sublasses.to.check){
        
        subclass.Data <- DataMatrix[subclasses.list[[subCL]], ]
        
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
    
    if(!include.RSD.matrix){
        RSD.results <- which(RSD.ok >= groups.ok.threshold)    
    } else{
        RSD.results <- list()
        RSD.results[[1]] <- which(RSD.ok >= groups.ok.threshold)
        RSD.results[[2]] <- RSD.matrix
    }
    
    
    return(RSD.results)
    
}