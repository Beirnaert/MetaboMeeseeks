#' parallel beakedpi peak picker
#'
#' function to parallelize the bakedpi algorithm
#'
#' @param DataMatrix A matrix with feature groups as columns and samples for rows. 
#' @param subclasses.list A list with either A) in each list entry a number of integers or 
#' numerics corresponding to the rows that form the subclasses to be checked, or B) in each 
#' list entry characters or factors corresponding to those found in the 'labels' attribute.
#' 
#' 
#' @return 
#' yammsList - alist with the data matrix, groupinfo and sample classes
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#'
#' @export
yamms.parallel = function(Datafiles, yamms.colData, nCPU = -1, mzRANGE, slicewidth = 100, sliceoverlap = 5, dbandwidth = c(0.01, 10), dgridstep = c(0.01, 1), outfileDens = NULL, dortalign = FALSE){
    
    

    
    
    
   Nslices = ceiling(diff(mzRANGE)/slicewidth)
   
   mz.seq.low  = as.integer(round(seq(from = mzRANGE[1]-sliceoverlap, to = mzRANGE[2]-sliceoverlap, length.out = 1+Nslices ) ) )
   mz.seq.high = as.integer(round(seq(from = mzRANGE[1]+sliceoverlap, to = mzRANGE[2]+sliceoverlap, length.out = 1+Nslices ) ) )
   mz.seq.low[1] = mzRANGE[1]
   mz.seq.high[1] = mzRANGE[1]
   mz.seq.low[length(mz.seq.low)] = mzRANGE[2]
   mz.seq.high[length(mz.seq.high)] = mzRANGE[2]
   
   if (nCPU == -1) {
       nCPU <- max(c(1,(parallel::detectCores(all.tests = FALSE, logical = TRUE) - 2)))
       nCPU <- min(Nslices, nCPU)
   }
   
   
   cl <- parallel::makeCluster(nCPU)
   doSNOW::registerDoSNOW(cl)
   performanceList <- list()
   parCounter <- NULL
   pb <- txtProgressBar(max=Nslices, style=3)
   progress <- function(n) setTxtProgressBar(pb, n)
   opts <- list(progress=progress)
   
   Slice.list <- foreach::foreach(parCounter = 1:Nslices, .options.snow =opts, .inorder = FALSE) %dopar% 
   {
       source("yamss/bakedpi.R")
       source("yamss/DataClasses.R")
       source("yamss/functions.R")
       source("yamss/reading.R")
       source("yamss/slicepi.R")
       source("yamss/utils.R")
     ptm = proc.time()  
    cmsRaw <- readMSdata(files = Datafiles, mzsubset = c(mz.seq.low[parCounter],mz.seq.high[parCounter+1]), colData = yamms.colData, verbose = FALSE)
       
    cmsProc <- bakedpi(cmsRaw, dbandwidth = c(0.01, 10), dgridstep = c(0.01, 1), mzsubset = c(mz.seq.low[parCounter],mz.seq.high[parCounter+1]), outfileDens = NULL, dortalign = FALSE,  verbose = FALSE)
    cmsSlice <- slicepi.light(cmsProc, cutoff = NULL, verbose = FALSE)
    
    cmsGroups = cmsSlice$rowData
    cmsMatrix = cmsSlice$assays$peakQuants 
    cmsClass = cmsSlice$colData$sampClasses
    proc.time() - ptm
    results = list(cmsGroups = cmsGroups, cmsMatrix = cmsMatrix, cmsClass = cmsClass)
    
    return(results)
   }
   close(pb)
   parallel::stopCluster(cl)
   
   
}