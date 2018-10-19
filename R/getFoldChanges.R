#' calculate fold changes
#'
#' Filter out isotopes based on the PUTMEDID concept
#'
#' @param DataMatrix The data matrix
#' @param classlabels the class labels (nr of vector elements must be equal to number of rows in DataMatrix)
#' @param baseClass The baseclass for the fold change clauculation
#' @param classOfInterest The class of interest for the fold change calculation
#' @param datatype Indication whether the Datamatrix is intensities (ints) or log2 values (log2s).
#'  
#' 
#'   
#' @return foldChanges A vector with the fold changes
#' 
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' 
#' 
#' @export
getFoldChanges <- function(DataMatrix, classlabels, baseClass = NULL, classOfInterest = NULL, datatype = "ints"){
    
    if(is.null(baseClass)){
        stop("No baseClass provided")
    }
    if(is.null(classOfInterest)){
        stop("No classOfInterest provided")
    }
    if(! baseClass %in% classlabels){
        stop("baseClass not in classlabels")
    }
    if(! classOfInterest %in% classlabels){
        stop("classOfInterest not in classlabels")
    }
    if(length(classlabels) != nrow(DataMatrix)){
        stop("Nr. of classlabels does not match the amount of rows in the DataMatrix")
    }
    
    if(datatype == "ints"){
        foldChanges = colMeans(DataMatrix[classlabels == classOfInterest,])/ colMeans(DataMatrix[classlabels == baseClass,])
    } else if(datatype == "log2s"){
        foldChanges = 2^(colMeans(DataMatrix[classlabels == classOfInterest,]) - colMeans(DataMatrix[classlabels == baseClass,]))
    } else{
        warning("no datatype provided, assumed data are intensities")
        foldChanges = colMeans(DataMatrix[classlabels == classOfInterest,])/ colMeans(DataMatrix[classlabels == baseClass,])
    }
    
    return(foldChanges)
}