#' NOT FINISHED!!!!!! Function to compare different RT corrections by looking at the most shifted peaks from in method
#'
#' Filter out isotopes based on the PUTMEDID concept
#'
#' @param RTcorr1.xcmsGrouped XCMS object of first RT correction method (post RT correction and subsequent grouping)
#' @param RTcorr2.xcmsGrouped XCMS object of second RT correction method (post RT correction and subsequent grouping)
#'  
#' 
#'   
#' @return 
#' 
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' 
#' 
#' @export
#' 
#' @importFrom purrr map2
#' 
#' 
RTcorrChecker = function(RTcorr1.xcmsGrouped = test.rtcorrs.grouped[[1]] , 
                         RTcorr2.xcmsGrouped = test.rtcorrs.obiwarp.grouped[[1]], 
                         plotTitle.method1 = "Loess.correction, span = 1", 
                         plotTitle.method2 = "Obiwarp correction, respons = 10",
                         RTdiff.cutoff.value = NULL,
                         RTdiff.cutoff.quantile = NULL){
    
    RT1.diffs = purrr::map2(.x = RTcorr1.xcmsGrouped@rt$corrected,
                            .y = RTcorr1.xcmsGrouped@rt$raw,
                            .f = ~ abs(.x - .y)
    )
    
    RT2.diffs = purrr::map2(.x = RTcorr2.xcmsGrouped@rt$corrected,
                            .y = RTcorr2.xcmsGrouped@rt$raw,
                            .f = ~ abs(.x - .y)
    )
    
    
    
    RTdiff.df = data.frame(RTdiffs = c(unlist(RT1.diffs), unlist(RT2.diffs)), 
                           RTmethod = c(rep("method1", length(unlist(RT1.diffs))),
                                        rep("method2", length(unlist(RT2.diffs)))) )
    
    ggplot(RTdiff.df, aes(RTdiffs,fill=RTmethod, colour = RTmethod))+ 
        geom_histogram(alpha = 0.5, position = "identity") + 
        theme_bw()
    
    
    if(is.null(RTdiff.cutoff) & is.null(RTdiff.cutoff.quantile)){
        warning("no RTdiff.cutoff or RTdiff.cutoff.quantile chosen, set to the 95 % quantile")
        RTdiff.cutoff = quantile(RTdiff.df$RTdiffs, probs = 0.95)
    } else if(is.null(RTdiff.cutoff)){
        if(RTdiff.cutoff.quantile > 1 | RTdiff.cutoff.quantile < 0){
            warning("RTdiff.cutoff.quantile not in [0,1] range. set to 0.95")
            RTdiff.cutoff.quantile = 0.95
        }
        RTdiff.cutoff = quantile(RTdiff.df$RTdiffs, probs = RTdiff.cutoff.quantile)
    }
    
    return(0)
    
    
}