#' Feature Quality score
#'
#' Calculate the feature quality score based on the presence of non failed gaussian fits.
#'
#'  
#' @param XCMSobject An xcmsSet object
#' @param indices.to.keep Vector with indices to keep. Either numeric indices or a vector with TRUE/FALSE (same length as nr. of groups in XCMSobject)
#' 
#'     
#' @return 
#' A vector with QC scores (between 0 and 1)
#' 
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#'
#'  
#' @export
Feature.quality.score = function(XCMSobject, indices.to.keep = NULL){

    if(! "xcmsSet" %in% class(XCMSobject)){
        stop("XCMSobject is not an xcmsSet")
    }
    
    if(length(XCMSobject@groupidx) == 0){
        stop("No groups in XCMSobject. Perform grouping step first")
    }
    
    if(is.null(indices.to.keep)){
        indices.to.keep = seq_along(XCMSobject@groupidx)
    } else if( "logical" %in% class(indices.to.keep)){
        indices.to.keep = which(indices.to.keep)
    }
    
    if(max(indices.to.keep) > length(XCMSobject@groupidx)){
        stop("The maximal index to keep is larger than the available amount of features.")
    }
    
    group.data = data.frame(XCMSobject@groups)
    
    
    Quality.score = vector("double",length(XCMSobject@groupidx))
    for(k in seq_along(indices.to.keep)){
        Quality.score = 
    }
    
    
    
}