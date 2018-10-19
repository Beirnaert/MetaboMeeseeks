#' Function to do peak filling (xcms) with restored RT values 
#'
#' 
#'
#' @param xcmsObject XCMS object to be peak filled
#' @param unStructureBatchesObject XCMS object after retreived after the unStructureBatches fumction
#'
#'   
#' @return xcmsObject A peak filled xcms object
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#' 
#' @importFrom xcms fillPeaks
#' 
#' @export
batchPeakFilling <- function(xcmsObject, unStructureBatchesObject){
    xcmsObject@rt$raw = unStructureBatchesObject@rt$Original_raw
    
    xcmsObject.peakfilled = xcms::fillPeaks(xcmsObject)
    
    return(xcmsObject.peakfilled)
}
