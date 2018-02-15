#' Function to filter out isotopes
#'
#' Filter out isotopes based on the PUTMEDID concept
#'
#' @param Data This can be either a) a matrix with 3 columns (mz, RT, Intensity) without missing values, b) an xcmsSet object post grouping (groups are features), or c) an xcmsSet@peaks object in which case each sample is filtered seperately (slowest)
#' @param max.ppm Maximal ppm distance
#' @param max.mz Maximal mz distance (not necessary if ppm supplied)
#' @param max.RT.narrow Prefered maximal RT distance 
#' @param max.RT.wide Extended RT search area. A peak is only accepted here when the improvement in ppm difference is at least accept.ppm.diff
#' @param accept.ppm.diff The improvement in ppm difference a peak in the wider RT area has to have for it to be accepted 
#' @param correlation Type of correlation method to be used. 
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
IsotopesBeGone = function(DataMatrix, subclasses.list = NULL, blanco.entry = NULL, blanco.data = NULL, Blanco.multiplier = 10,ignore.missing.values = TRUE, labels = NULL, feature.orientation = "columns", groups.ok.threshold = 1){
    
    
    
    
}