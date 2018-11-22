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
#' @param max.RT.width The maximal allowed RT width for a group
#' @param max.ppm.diff The maximal allowed ppm difference within a group
#' @param max.dupl.prop The maximal allowed proportion of duplicated peaks in the group
#' @param max.missing.prop The maximal allowed proportion of missing peaks in the group
#' @param thresholds.to.pass The number of thresholds to pass (default is all)
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
group.filter= function(XCMSobject, max.RT.width = NULL, max.ppm.diff = NULL, max.dupl.prop = NULL, max.missing.prop = NULL, thresholds.to.pass = NULL){
    
    if(is.null(thresholds.to.pass)){
        thresholds.to.pass = sum(!c(is.null(max.RT.width), is.null(max.ppm.diff), is.null(max.dupl.prop), is.null(max.missing.prop)))
    }
    
    #if( WithinGroupsOfInterest & is.null(ClassesOfInterest)){
    #    stop("ClassesOfInterest is not provided. This is necessary when WithinGroupsOfInterest = TRUE")
    #}
    
    #if( !is.null(ClassesOfInterest) & is.null(classVector)){
    #    stop("classVector is not provided. This is necessary when ClassesOfInterest is provided")
    #}
    
    
    
    nGroups = length(XCMSobject@groupidx)
    group.stats = matrix(NA,ncol = 8, nrow = nGroups)
    group.stats = data.frame(group.stats)
    colnames(group.stats) = c("mzmed", "RTmed", "ppm.diff", "RT.diff", "prop.missing", "nonMissingPeaks", "nMissing", "prop.dupl")
    
    
    groups.df = data.frame(XCMSobject@groups)
    
    group.stats$mzmed = groups.df$mzmed
    group.stats$RTmed = groups.df$rtmed
    
    group.stats$ppm.diff = 10^6 * (groups.df$mzmax - groups.df$mzmin) / groups.df$mzmed
    group.stats$RT.diff =  groups.df$rtmax - groups.df$rtmin
    
    group.stats$nonMissingPeaks = groups.df$npeaks
    
    for(k in 1:nGroups){
        
        groupPeaks = data.frame(XCMSobject@peaks[XCMSobject@groupidx[[k]],])
      
        group.stats$prop.missing[k] = sum(is.na(groupPeaks$intb)) / nrow(groupPeaks)
        
        group.stats$nMissing[k] = sum(is.na(groupPeaks$intb))
        
        group.stats$prop.dupl[k] = sum(duplicated(groupPeaks$sample)) / groups.df$npeaks[k]
        
    }
    
    group.stats$score = 0
    if(!is.null(max.RT.width)){
        group.stats$score[group.stats$RT.diff <= max.RT.width] = group.stats$score[group.stats$RT.diff <= max.RT.width] + 1
    }
    if(!is.null(max.RT.width)){
        group.stats$score[group.stats$ppm.diff <= max.ppm.diff] = group.stats$score[group.stats$ppm.diff <= max.ppm.diff] + 1
    }
    if(!is.null(max.RT.width)){
        group.stats$score[group.stats$prop.dupl <= max.dupl.prop] = group.stats$score[group.stats$prop.dupl <= max.dupl.prop] + 1
    }
    if(!is.null(max.RT.width)){
        group.stats$score[group.stats$prop.missing <= max.missing.prop] = group.stats$score[group.stats$prop.missing <= max.missing.prop] + 1
    }
    
    
    XCMSobject.filtered = XCMSobject
    XCMSobject.filtered@groups = XCMSobject@groups[group.stats$score >= thresholds.to.pass, ]
    XCMSobject.filtered@groupidx = XCMSobject@groupidx[group.stats$score >= thresholds.to.pass]
    
    return(XCMSobject.filtered)
    
}