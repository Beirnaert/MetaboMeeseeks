#' CAMERA.accumulator
#'
#' Function to accumulate CAMERA results for all combinations of given parameters.
#'
#'  
#' @param XCMS.object xcmsSet object
#' @param polarity xcmsSet polarity
#' @param perfwhm xcms perfwhm parameter
#' @param ppm_threshold xcms ppm_threshold parameter
#' @param ppm_threshold.multiplier parameter with which the ppm_threshold is multiplied
#' @param mzabs xcms mzabs parameter
#' @param intval xcms intval parameter
#' @param maxcharge xcms maxcharge parameter
#' @param nCPU The number of CPUs to use. Default is maximum available minus 2.
#' @param Standards.list Yet to implement
#' 
#' @return 
#'    
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' 
#' @importFrom foreach %dopar% foreach
#' 
#' @export
#' 
CAMERA.accumulator = function(XCMS.object, polarity = NULL, perfwhm = 0.6, ppm_threshold = 10, ppm_threshold.multiplier = 1.5, mzabs = 0.02, intval = "into", maxcharge = 2, nCPU = -1, Standards.list=NULL){   
    
    
    if(!polarity %in% c("negative", "positive") ){
        stop("The 'polarity' parameter is not one of 'negative' or 'positive")
    }
    
    Combos = expand.grid(perfwhm = perfwhm, 
                         ppm_threshold = ppm_threshold, 
                         ppm_threshold.multiplier = ppm_threshold.multiplier,
                         mzabs = mzabs,
                         intval = intval,
                         maxcharge = maxcharge, 
                         stringsAsFactors = FALSE)
    
    nCombos = nrow(Combos)
    
    if (nCPU == -1) {
        nCPU <- max(c(1,(parallel::detectCores(all.tests = FALSE, logical = TRUE) - 2)))
        nCPU <- min(nCombos, nCPU)
    }
    
    
    cl <- parallel::makeCluster(nCPU)
    doSNOW::registerDoSNOW(cl)
    performanceList <- list()
    parCounter <- NULL
    print("Running CAMERA iterations.")
    pb <- txtProgressBar(max=nCombos, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    
    # divide the number of runs for each core
    
    CAMERA.iters <- foreach::foreach(parCounter = 1:nCombos, .options.snow=opts, .inorder = TRUE, .packages = c("CAMERA")) %dopar%
    {
        
        an_neg <- xsAnnotate(XCMS.object, polarity = polarity)
        an_neg2 <- groupFWHM(an_neg, perfwhm = Combos$perfwhm[parCounter]) # group peaks by retention time
        an_neg3 <- findIsotopesWithValidation(object = an_neg2, ppm = Combos$ppm_threshold[parCounter] * Combos$ppm_threshold.multiplier[parCounter], mzabs = Combos$mzabs[parCounter], intval = Combos$intval[parCounter] , maxcharge = Combos$maxcharge[parCounter]) # annotate isotopic peaks
        
        Npseudospec <- length(an_neg2@pspectra)
        Nisotopes <- nrow(an_neg3@isoID)
        
        Present = data.frame(matrix(0,nrow = 1, ncol = nrow(XCMS.object@groups)))
        colnames(Present) = paste("group_",as.character(seq(1,nrow(XCMS.object@groups))))
        Present[,an_neg3@isoID[,1]] = an_neg3@isoID[,3]
        
        Results = data.frame(Combos[parCounter,], Npseudospec, Nisotopes, Present)
        
        return(Results)
    }
    close(pb)
    parallel::stopCluster(cl)
    
    CAMERA.iters <- data.table::rbindlist(CAMERA.iters)
    
    
    overview.plot <- ggplot() +
                        geom_line(data = CAMERA.iters , aes(x=perfwhm, y = Npseudospec, colour = "pseudospec")) +
                        geom_line(data = CAMERA.iters , aes(x=perfwhm, y = Nisotopes, colour = "isotopes")) +
                        theme_bw()
    
    print(overview.plot)
    
    
    
    # the standards matching
    
    for(std in 1:length(Standards.list)){
        
    }
    
    
    
    
    # The image plot
    iso.groups <- as.matrix(CAMERA.iters[,9:ncol(CAMERA.iters)])
    
    nLevels <- max(iso.groups)
    angles.col <- seq(0, 360, length.out = nLevels + 1) + 15
    plotcols <- grDevices::hcl(h = angles.col, l = 65, c = 120)[1:nLevels]
    
    
    image(iso.groups, zlim = c(0.5,max(iso.groups)+1), col = plotcols, breaks = 0.5+seq(0,ceiling(max(iso.groups))))
    
    # The profile plot
    
    # the ones that are constant are okay
    non.changing.iso <- matrixStats::colMins(iso.groups) == matrixStats::colMaxs(iso.groups) & matrixStats::colMaxs(iso.groups) > 0
    non.changing.0 <- matrixStats::colMins(iso.groups) == matrixStats::colMaxs(iso.groups) & matrixStats::colMaxs(iso.groups) == 0
    changing.from0 <- matrixStats::colMins(iso.groups) != matrixStats::colMaxs(iso.groups) & matrixStats::colMins(iso.groups) == 0
    changing.from_non0 <- matrixStats::colMins(iso.groups) != matrixStats::colMaxs(iso.groups) & matrixStats::colMins(iso.groups) > 0
    

    
    profiles.df.non.changing.iso = data.frame(X = rep(CAMERA.iters$perfwhm, ncol(iso.groups[,non.changing.iso]) ), Y= as.numeric(iso.groups[,non.changing.iso]), group = as.numeric(matrix(rep(seq(1,ncol(iso.groups[,non.changing.iso])),nrow(CAMERA.iters)), ncol = ncol(iso.groups[,non.changing.iso]), nrow = nrow(CAMERA.iters), byrow = T)))
    profiles.df.changing.from0 = data.frame(X = rep(CAMERA.iters$perfwhm, ncol(iso.groups[,changing.from0]) ), Y= as.numeric(iso.groups[,changing.from0]), group = as.numeric(matrix(rep(seq(1,ncol(iso.groups[,changing.from0])),nrow(CAMERA.iters)), ncol = ncol(iso.groups[,changing.from0]), nrow = nrow(CAMERA.iters), byrow = T)))
    profiles.df.changing.from_non0 = data.frame(X = rep(CAMERA.iters$perfwhm, ncol(iso.groups[,changing.from_non0]) ), Y= as.numeric(iso.groups[,changing.from_non0]), group = as.numeric(matrix(rep(seq(1,ncol(iso.groups[,changing.from_non0])),nrow(CAMERA.iters)), ncol = ncol(iso.groups[,changing.from_non0]), nrow = nrow(CAMERA.iters), byrow = T)))
    
    ggplot(profiles.df.non.changing.iso, aes(x=X,y=Y, group = group, colour =group)) +
        geom_line() + 
        theme_bw () + 
        ggtitle(paste("nr. of isotopoplogues = ",as.character(max(profiles.df.non.changing.iso$group)))) 
    
    ggplot(profiles.df.changing.from0, aes(x=X,y=Y, group = group, colour =group)) +geom_line() +
        geom_line() + 
        theme_bw () + 
        ggtitle(paste("nr. of isotopoplogues = ",as.character(max(profiles.df.changing.from0$group)))) 
    
    ggplot(profiles.df.changing.from_non0, aes(x=X,y=Y, group = group, colour =group)) +geom_line() +
        geom_line() + 
        theme_bw () + 
        ggtitle(paste("nr. of isotopoplogues = ",as.character(max(profiles.df.changing.from_non0$group)))) 
}