#' Group Filtering
#'
#' After grouping and peak filling, this function filters out bad groups/features. 
#' There is the possibility of supplying groups of interest in which case the. 
#' TODO: param ToDo......ClassesOfInterest (optional) The classes to include in the calculation 
#' of group properties. 
#' param ToDo classVector (required if ClassesOfInterest provided) vector with classlabels.
#' param ToDo.......WithinGroupsOfInterest (default is FALSE), Option to choose to do the 
#' calculations within class instead of over the entire group/feature
#'
#'  
#' @param XCMSobject An xcmsSet object after peak picking.
#' @param max_RT_width The maximal allowed RT width for a group
#' @param max_ppm_diff The maximal allowed ppm difference within a group
#' @param max_dupl_prop The maximal allowed proportion of duplicated peaks in the group
#' @param max_missing_prop The maximal allowed proportion of missing peaks in the group
#' @param thresholds_to_pass The number of thresholds to pass (default is all)
#'  
#' 
#' 
#' 
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#'
#'  
#' @export
GroupFiltR= function(XCMSobject, max_RT_width = NULL, max_ppm_diff = NULL, max_dupl_prop = NULL, max_missing_prop = NULL, thresholds_to_pass = NULL){
    
    if(is.null(thresholds_to_pass)){
        thresholds_to_pass <- sum(!c(is.null(max_RT_width), is.null(max_ppm_diff), is.null(max_dupl_prop), is.null(max_missing_prop)))
    }
    
    
    nGroups <- length(XCMSobject@groupidx)
    group_stats <- matrix(NA,ncol = 8, nrow = nGroups)
    group_stats <- data.frame(group_stats)
    colnames(group_stats) <- c("mzmed", "RTmed", "ppm_diff", "RT_diff", "prop_missing", "nonMissingPeaks", "nMissing", "prop_dupl")
    
    
    groups_df <- data.frame(XCMSobject@groups)
    
    group_stats$mzmed <- groups_df$mzmed
    group_stats$RTmed <- groups_df$rtmed
    
    group_stats$ppm_diff <- 10^6 * (groups_df$mzmax - groups_df$mzmin) / groups_df$mzmed
    group_stats$RT_diff <-  groups_df$rtmax - groups_df$rtmin
    
    group_stats$nonMissingPeaks <- groups_df$npeaks
    
    for(k in 1:nGroups){
        
        groupPeaks <- data.frame(XCMSobject@peaks[XCMSobject@groupidx[[k]],])
      
        group_stats$prop_missing[k] <- sum(is.na(groupPeaks$intb)) / nrow(groupPeaks)
        
        group_stats$nMissing[k] <- sum(is.na(groupPeaks$intb))
        
        group_stats$prop_dupl[k] <- sum(duplicated(groupPeaks$sample)) / groups_df$npeaks[k]
        
    }
    
    group_stats$score <- 0
    if(!is.null(max_RT_width)){
        group_stats$score[group_stats$RT_diff <= max_RT_width] <- group_stats$score[group_stats$RT_diff <= max_RT_width] + 1
    }
    if(!is.null(max_RT_width)){
        group_stats$score[group_stats$ppm_diff <= max_ppm_diff] <- group_stats$score[group_stats$ppm_diff <= max_ppm_diff] + 1
    }
    if(!is.null(max_RT_width)){
        group_stats$score[group_stats$prop_dupl <= max_dupl_prop] <- group_stats$score[group_stats$prop_dupl <= max_dupl_prop] + 1
    }
    if(!is.null(max_RT_width)){
        group_stats$score[group_stats$prop_missing <= max_missing_prop] <- group_stats$score[group_stats$prop_missing <= max_missing_prop] + 1
    }
    
    
    XCMSobject_filtered <- XCMSobject
    XCMSobject_filtered@groups <- XCMSobject@groups[group_stats$score >= thresholds_to_pass, ]
    XCMSobject_filtered@groupidx <- XCMSobject@groupidx[group_stats$score >= thresholds_to_pass]
    
    return(XCMSobject_filtered)
    
}