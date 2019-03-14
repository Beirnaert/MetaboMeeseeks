#' Standards Locater
#'
#' This function searches for XCMS objects (peaks, groups) that match with the provided standards (optionally with isotopic profiles).
#'
#'  
#' @param Standards.mz A vector (single standards) or matrix (isotopic profiles) with the standards.
#' @param Standards.RT A vector with the RT values
#' @param Data This can be either a) a matrix with 2 or 3* columns (mz, RT, Intensity*) without missing values, b) an xcmsSet object post grouping (the groups are features, intensity ratio's can be checked) or c) an XCMS@groups matrix object (if only this is supplied there will be no intensity ratio check). Note that an XCMS@peaks object is also supported but the function only returns 1 match per standard (the closest one, or the first of a number of equally close ones)
#' @param Standards.ratios Optional matrix. Matrix with the intensity ratio's
#' @param max.ppm maximal ppm distance
#' @param max.mz maximal mz distance (not necessary if ppm supplied)
#' @param max.RT.narrow preferred maximal RT distance 
#' @param max.RT.wide extended RT search area. A peak is only accepted here when the improvement in ppm difference is at least accept.ppm.diff
#' @param accept.ppm.diff The improvement in ppm difference a peak in the wider RT area has to have for it to be accepted 
#'  
#' @return 
#' A list with 3 vector/matrix entries (the same size as the supplied vector/matrix with standards). 
#' The first are the matched peaks/groups. The second the ppm differences. The third the RT differences.
#' 
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#'
#' @export
StdMatchR = function(Standards.mz, Standards.RT, Data, Standards.ratios = NULL, max.ppm = NULL, max.mz = NULL, max.RT.narrow = NULL, max.RT.wide = NULL, accept.ppm.diff = NULL){
    
    ratioCheck = FALSE
    
    if( !is.null(Standards.ratios) ){
        if( ! "matrix" %in% class(Standards.ratios) ){
            stop("Standards.ratios is not a matrix. Should be of format (Standard X Intensitie ratio)")
        }
        if(!"xcmsSet" %in% class(Data) ){
            if(ncol(Data) != 3){
                stop("Data is not an xcmsSet object, nor is it a matrix with 3 columns (mz, RT, Intensity). One of these is necassary for the ratio check")
            }
        }
        if( !all( dim(Standards.ratios) ==  dim(Standards.mz)) ) {
            stop("The dimensions of Standards.ratios is not equal to those of Standards.mz")
        }
        ratioCheck = TRUE
    }
    
    if("matrix" %in% class(Standards.mz)){
        n.stds = dim(Standards.mz)[1]  
        n.isotopes = dim(Standards.mz)[2]    
    } else if("numeric" %in% class(Standards.mz)){
        n.stds = length(Standards.mz)  
        n.isotopes = 1
        Standards.mz = matrix(Standards.mz, nrow = n.stds, ncol = n.isotopes)
    }
    
    xcmsRatio = FALSE
    if("xcmsSet" %in% class(Data)){
        xcmsData = Data
        Data = xcmsData@groups[,c(1,4)]
        xcmsRatio = TRUE
        
    } else if(ncol(Data) > 3){
        Data = Data[,c(1,4)]
        ratioCheck = FALSE
    } else if(ncol(Data) == 2){
        ratioCheck = FALSE
    }
    
    if( is.null(max.ppm) & is.null(max.mz) ){
        warning("No maximal values submitted for the mz or ppm difference. max.ppm set to 50 and max.mz to 0.2")
        max.ppm = 50
        max.mz = 0.2
    }
    
    if(is.null(max.mz)){
        max.mz = max(Data[,1])*max.ppm/(10^6) + 0.1
    }
    
    if(is.null(max.ppm)){
        max.ppm = 10^6 * max.mz/min(Data[,1])
    }
    
    if(is.null(max.RT.narrow)){
        max.RT.narrow = 300 # 5 minutes
    }
    if(is.null(max.RT.wide)){
        max.RT.wide = max.RT.narrow 
    }
    
    if(is.null(accept.ppm.diff)){
        accept.ppm.diff = 5
    }
    
    standards_groups = matrix(NA, nrow = n.stds, ncol = n.isotopes)
    standards_groups.onlynarrow = matrix(NA, nrow = n.stds, ncol = n.isotopes)
    standards_groups.ppmDiff = matrix(NA, nrow = n.stds, ncol = n.isotopes)
    standards_groups.RTDiff = matrix(NA, nrow = n.stds, ncol = n.isotopes)
    
    for(mt in 1:n.stds){
        
        for(is in 1:n.isotopes){
            
            if( !is.na(Standards.mz[mt,is])){
                mz.match = abs(Data[,1] - Standards.mz[mt,is])
                ppm.match = 10^6 * abs(Data[,1] - Standards.mz[mt,is])/Standards.mz[mt,is]
                RT.match = abs(Data[,2] - Standards.RT[mt])
                
                matched.mz = mz.match < max.mz
                matched.ppm = ppm.match < max.ppm
                matched.RT.narrow = RT.match < max.RT.narrow
                matched.RT.wide = RT.match < max.RT.wide
                
                if(sum(matched.mz & matched.ppm & matched.RT.narrow ) > 0){
                    
                    ranking.ppmfirst.narrow = order(ppm.match[matched.mz & matched.ppm & matched.RT.narrow], RT.match[matched.mz & matched.ppm & matched.RT.narrow], decreasing = FALSE)[1]
                    ranking.ppmfirst.wide = order(ppm.match[matched.mz & matched.ppm & matched.RT.wide], RT.match[matched.mz & matched.ppm & matched.RT.wide], decreasing = FALSE)[1]
                    
                    top.peakindex.narrow = which( matched.mz & matched.ppm & matched.RT.narrow >0)[ranking.ppmfirst.narrow]
                    top.peakindex.wide = which( matched.mz & matched.ppm & matched.RT.wide >0)[ranking.ppmfirst.wide]
                    
                    standards_groups.onlynarrow[mt,is] = top.peakindex.narrow
                    
                    if( ppm.match[top.peakindex.wide] < (ppm.match[top.peakindex.narrow] - accept.ppm.diff)  ){
                        standards_groups[mt,is] = top.peakindex.wide
                        standards_groups.ppmDiff[mt,is] = ppm.match[top.peakindex.wide]
                        standards_groups.RTDiff[mt,is] = RT.match[top.peakindex.wide]
                    } else{
                        standards_groups[mt,is] = top.peakindex.narrow
                        standards_groups.ppmDiff[mt,is] = ppm.match[top.peakindex.narrow]
                        standards_groups.RTDiff[mt,is] = RT.match[top.peakindex.narrow]
                    }
                    
                    #if(length(closest) > 1){
                    #    closets.
                    #}
                    #standards_groups[mt,is] = which( sum(matched.mz & matched.ppm & matched.RT) >0)[ranking.ppmfirst]
                } else if(sum(matched.mz & matched.ppm & matched.RT.wide ) > 0){
                    ranking.ppmfirst.wide = order(ppm.match[matched.mz & matched.ppm & matched.RT.wide], RT.match[matched.mz & matched.ppm & matched.RT.wide], decreasing = FALSE)[1]
                    top.peakindex.wide = which( matched.mz & matched.ppm & matched.RT.wide >0)[ranking.ppmfirst.wide]
                    
                    standards_groups[mt,is] = top.peakindex.wide
                    standards_groups.ppmDiff[mt,is] = ppm.match[top.peakindex.wide]
                    standards_groups.RTDiff[mt,is] = RT.match[top.peakindex.wide]
                }
            }
        }
        
    }
    
    if(ratioCheck){
        if(!xcmsRatio){
            for(mt2 in 1:n.stds){
                present = !is.na(standards_groups[mt2,])
                present.intensities = Data[standards_groups[mt2,present],3]
                present.ratios = Standards.ratios[mt2,present]
                
                passed.lower = rep(TRUE, length(present.intensities))
                for(jj in 2:length(present.intensities)){
                    if(present.intensities[jj] >= min(present.intensities[1:(jj-1)]) ){
                        passed.lower[jj] = FALSE
                    }
                }
                
                standards_groups[mt2,present][!passed.lower] = NA
                standards_groups.ppmDiff[mt2,present][!passed.lower] = NA
                standards_groups.RTDiff[mt2,present][!passed.lower] = NA
            }
        } else if (xcmsRatio){
            plot.output = list()
            for(mt2 in 1:n.stds){
                
                # here we calculate the ratio based on 1 sample, appearing in all features.
                present = !is.na(standards_groups[mt2,])
                present.index = standards_groups[mt2,present]
                present.ratios = Standards.ratios[mt2,present]
                
                samples.in.feature = list()
                for(smp in 1:length(present.index)){
                    sample.peakdata = data.frame(xcmsData@peaks[xcmsData@groupidx[[present.index[smp] ]], ])
                    samples.in.feature[[smp]] = sample.peakdata$sample[!is.na(sample.peakdata$intb)]
                }
                
                possible.samples = Reduce(intersect, samples.in.feature)
                if(length(possible.samples) == 0){
                    possible.samples.sub = rep(NA,(length(present.index)-1))
                    for(red in 1:(length(present.index)-1) )
                        possible.samples.sub[red] = length(Reduce(intersect, samples.in.feature[1:(length(present.index) - red)]))
                    
                    
                    to.reduce = which(possible.samples.sub>0)
                    
                    
                    standards_groups[mt2,present][(length(present.index) - to.reduce + 1):length(present.index)] = NA
                    standards_groups.ppmDiff[mt2,present][(length(present.index) - to.reduce + 1):length(present.index)] = NA
                    standards_groups.RTDiff[mt2,present][(length(present.index) - to.reduce + 1):length(present.index)] = NA
                    
                    present = !is.na(standards_groups[mt2,])
                    present.index = standards_groups[mt2,present]
                    present.ratios = Standards.ratios[mt2,present]
                    
                    samples.in.feature = list()
                    for(smp in 1:length(present.index)){
                        sample.peakdata = data.frame(xcmsData@peaks[xcmsData@groupidx[[present.index[smp] ]], ])
                        samples.in.feature[[smp]] = sample.peakdata$sample[!is.na(sample.peakdata$intb)]
                    }
                    
                }
                
                
                
                if(length(possible.samples)!=0){
                    
                    possible.samples.peakdata = data.frame(xcmsData@peaks[xcmsData@groupidx[[present.index[1] ]], ])
                    possible.samples.peakdata = possible.samples.peakdata[possible.samples.peakdata$sample %in% possible.samples,c("into","sample")]
                    possible.samples.order = order(possible.samples.peakdata$into, decreasing = TRUE)
                    possible.samples.order = possible.samples.order[!duplicated(possible.samples.peakdata$sample[possible.samples.order])]
                    
                    start = 1
                    stop = 0
                    all.ratio.int = matrix( NA, ncol = 4, nrow = length(possible.samples)*length(present.index))
                    for(jj in 1:length(possible.samples.order)){
                        
                        chosen.sample = possible.samples.peakdata$sample[possible.samples.order[jj]]
                        
                        present.intensities = rep(NA, length(present.index))
                        for(smp in 1:length(present.index)){
                            sample.peakdata = data.frame(xcmsData@peaks[xcmsData@groupidx[[present.index[smp] ]], ])
                            present.intensities[smp] = max(sample.peakdata$into[sample.peakdata$sample == chosen.sample])[1]
                        }
                        
                        passed.lower = rep(TRUE, length(present.intensities))
                        if(length(present.intensities) > 1){
                            for(jj in 2:length(present.intensities)){
                                if(present.intensities[jj] >= min(present.intensities[1:(jj-1)]) ){
                                    passed.lower[jj] = FALSE
                                }
                            }
                        }
                        
                        standards_groups[mt2,present][!passed.lower] = NA
                        standards_groups.ppmDiff[mt2,present][!passed.lower] = NA
                        standards_groups.RTDiff[mt2,present][!passed.lower] = NA
                        
                        stop=stop+length(present.index)
                        all.ratio.int[start:stop,1] = present.ratios
                        all.ratio.int[start:stop,2] = present.intensities
                        all.ratio.int[start:stop,3] = chosen.sample
                        all.ratio.int[start:stop,4] = as.numeric(passed.lower)
                        start = stop + 1
                    }
                    
                    
                    
                    all.ratio.int = data.frame(all.ratio.int)
                    colnames(all.ratio.int) = c("ratio", "intensity", "sample", "IR.check")
                    all.ratio.int$sample = as.factor(all.ratio.int$sample)
                    
                    plot.output[[mt2]] <- ggplot() + 
                        geom_line(data = subset(all.ratio.int, ~IR.check ==1), 
                                  aes(x = ~ratio, y = ~intensity, colour = ~sample)) + 
                        geom_point(data = subset(all.ratio.int, ~IR.check ==1), 
                                   aes(x = ~ratio, y = ~intensity, colour = ~sample))+
                        geom_point(data = subset(all.ratio.int, ~IR.check ==0), 
                                   aes(x = ~ratio, y = ~intensity, colour = ~sample), shape = 4)+
                        scale_x_reverse()+
                        theme_bw() +
                        ggtitle(paste("Standard ",as.character(mt2))) +
                        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
                    
                }
                
                
            }
        }
        
    }
    
    
    #for(mt in 1:n.stds){
    #    first.na = which(is.na(standards_groups[mt,]))[1]
    #    standards_groups[mt,first.na:n.isotopes] = NA
    #}
    if(xcmsRatio & ratioCheck){
        Results = list(Matched = standards_groups, ppm.diffs = standards_groups.ppmDiff, RT.diffs = standards_groups.RTDiff, Plots = plot.output)
    } else{
        Results = list(Matched = standards_groups, ppm.diffs = standards_groups.ppmDiff, RT.diffs = standards_groups.RTDiff)
    }
    return(Results)
}