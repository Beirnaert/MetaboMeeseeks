#' Peak Filtering
#'
#' After peak picking this function filters out bad groups/features.
#'
#'  
#' @param XCMSobject An xcmsSet object after peak picking.
#' @param min.RT.diff The minimal retentiuon time width a peak must have
#' @param max.RT.diff The maximal retention time width a peak may have
#' @param max.ppm.diff The maximal allowed ppm width for a single peak
#' @param gauss.fail.filter Whether to filter out 
#' 
#' @return filtered xcmsSet object
#' 
#' 
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#'

#'  
#' @export
peak.filter= function(XCMSobject, min.RT.width = NULL, max.RT.width = NULL, max.ppm.width = NULL, gauss.fail.threshold = NULL){
 
    nPeaks = nrow(XCMSobject@peaks)
    peaks.df = data.frame(XCMSobject@peaks)
    
    RT.width = peaks.df$rtmax - peaks.df$rtmin
    ppm.width = 10^6 * (peaks.df$mzmax - peaks.df$mzmin) / peaks.df$mz
    
    pass.indices = rep(TRUE, nPeaks)
    
    if(!is.null(min.RT.width)){
        pass.indices = pass.indices & RT.width >= min.RT.width
    }
    
    if(!is.null(max.RT.width)){
        pass.indices = pass.indices & RT.width <= max.RT.width
    }
    
    if(!is.null(max.ppm.width)){
        pass.indices = pass.indices & ppm.width <= max.ppm.width
    }

    if(!is.null(gauss.fail.threshold)){
        pass.indices = pass.indices & !is.na(peaks.df$egauss) 
        pass.indices = pass.indices & peaks.df$egauss <= gauss.fail.threshold
    }
   
   
    XCMSobject.filtered = XCMSobject
    XCMSobject.filtered@peaks = XCMSobject@peaks[pass.indices,]
    
    return(XCMSobject.filtered)
    
}