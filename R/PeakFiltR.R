#' Peak Filtering
#'
#' After peak picking this function filters out bad groups/features.
#'
#'  
#' @param XCMSobject An xcmsSet object after peak picking.
#' @param min_RT_width The minimal retention time width a peak must have
#' @param max_RT_width The maximal retention time width a peak may have
#' @param max_ppm_width The maximal allowed ppm width for a single peak
#' @param gauss_fail_threshold The threshold for the egauss XCMS parameter, below which the peak is rejected.
#' 
#' @return filtered xcmsSet object
#' 
#' 
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#'  
#' @export
PeakFiltR <- function(XCMSobject, min_RT_width = NULL, max_RT_width = NULL, max_ppm_width = NULL, gauss_fail_threshold = NULL){
 
    nPeaks <- nrow(XCMSobject@peaks)
    peaks_df <- data.frame(XCMSobject@peaks)
    
    RT_width <- peaks_df$rtmax - peaks_df$rtmin
    ppm_width <- 10^6 * (peaks_df$mzmax - peaks_df$mzmin) / peaks_df$mz
    
    pass_indices <- rep(TRUE, nPeaks)
    
    if(!is.null(min_RT_width)){
        pass_indices <- pass_indices & RT_width >= min_RT_width
    }
    
    if(!is.null(max_RT_width)){
        pass_indices <- pass_indices & RT_width <= max_RT_width
    }
    
    if(!is.null(max_ppm_width)){
        pass_indices <- pass_indices & ppm_width <= max_ppm_width
    }

    if(!is.null(gauss_fail_threshold)){
        pass_indices <- pass_indices & !is.na(peaks_df$egauss) 
        pass_indices <- pass_indices & peaks_df$egauss <= gauss_fail_threshold
    }
   
   
    XCMSobject_filtered <- XCMSobject
    XCMSobject_filtered@peaks <- XCMSobject@peaks[pass_indices,]
    
    return(XCMSobject_filtered)
    
}