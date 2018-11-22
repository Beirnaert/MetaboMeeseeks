#' Filling and grouping checker
#'
#' Function to check the peak filling and grouping. Every group gets a score for the mz range, the RT range and the filling/
#' For example, filled peaks cannot have intensities considerably larger than the already present peaks. 
#' This violates the assumption that these peaks are missed in peak detection because they are smaller.
#'
#'  
#' @param XCMSobject An xcmsSet after peak grouping and optionally after peak filling.
#'   
#' @return 
#' A data.frame with quality control metrics for each group.        
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#'
#' @export
GroupFillCheck = function(XCMSobject){
    
    # check for nr duplicate samples,
    # ratio of peak filled samples 
    # distribution of peak filled peaks vs present peaks
    # spread in mz (ppm)
    # width in RT
    
    GroupFill.QC = data.frame(matrix(NA, ncol = 8))
    
    for(ft in 1:nrow(XCMSobject@groups)){
        grp = XCMSobject
    }
    
    
    
}