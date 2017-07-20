#' calculate fold changes
#'
#' Filter out isotopes based on the PUTMEDID concept
#'
#' @param RTcorr1.xcmsGrouped XCMS object of first RT correction method (post RT correction and subsequent grouping)
#' @param RTcorr2.xcmsGrouped XCMS object of second RT correction method (post RT correction and subsequent grouping)
#'  
#' 
#'   
#' @return foldChanges A vector with the fold changes
#' 
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' 
#' 
#' @export
getFoldChanges <- function(DataMatrix, classlabels, baseClass = NULL, classOfInterest = NULL){
    
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
    
    foldChanges = (colMeans(DataMatrix[classlabels == classOfInterest,]) - colMeans(DataMatrix[classlabels == baseClass,]))/colMeans(DataMatrix[classlabels == baseClass,])
    
    return(foldChanges)
}