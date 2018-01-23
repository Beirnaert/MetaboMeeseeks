#' Function to split up xcms object by catch to perform a within batch analysis and merge the results at the end
#'
#' binnen batch analyse 
#' gebruik obiwarp om te weten hoe je de features van batch 1 en 2 aan elkaar moet plakken
#' doe dan analyse loess, PLSDA van batch 1 en 2 afzonderlijk. Giet alles in data frame (significant, fold change)
# 
#'
#' @param xcmsObject XCMS object after multi-batch retention time correction and grouping.
#' @param BatchClassLabel The name of the vector containing the Batch identifier (in xcmsObject@phenoData)
#' @param BatchNames (optional) A vector with the names of the individual batches to consider eg c("Batch1", "Batch2"). Note that this is only if a subset of all the batches have to be selected. Otherwise supplying BatchClassLabel is sufficient. 
#'
#'   
#' @return xcmsObject.split A list with in each list entry an xcms object of a single batch.
#' 
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @importFrom utils head
#'  
#' 
#' @export
BatchSplitR <- function(xcmsObject, BatchClassLabel = NULL, BatchNames = NULL ){
    
    if(!"xcmsSet" %in% class(xcmsObject)){
        stop("xcmsObject is not of class 'xcmsSet'")
    }
    
    if(is.null(BatchClassLabel)){
        print(head(xcmsObject@phenoData))
        BatchClassLabel <- readline(paste("Which of the following classes describes the batches: ",paste(colnames(xcmsObject@phenoData), collapse=" or ")))
    }
    
    if(is.null(BatchNames)){
        BatchNames = as.character(unique(xcmsObject@phenoData[BatchClassLabel])[,1])
    } else if(!all(BatchNames %in% as.character(unique(xcmsObject@phenoData[BatchClassLabel])[,1]))) {
        stop(paste("The supplied BatchNames could not all be matched to items in xcmsObject@phenoData which consists of:", paste(as.character(unique(xcmsObject@phenoData[BatchClassLabel])[,1]), collapse=" and ")))
    }
    
    xcmsObject.split = list()
    #mutual.RT.correction = list()
    for(k in 1:length(BatchNames)){
        split.samples = which(xcmsObject@phenoData[BatchClassLabel] == BatchNames[k])
        
        splitObject = xcmsObject
        sampleColumn = which(colnames(xcmsObject@peaks) == "sample")
        splitObject@peaks = splitObject@peaks[as.data.frame(xcmsObject@peaks)$sample %in% split.samples, ]
        splitObject@peaks[,sampleColumn] = splitObject@peaks[,sampleColumn] - split.samples[1] + 1
        splitObject@groups = matrix(NA,ncol=0,nrow=0)
        splitObject@groupidx = list()
        phenoDat = xcmsObject@phenoData[split.samples,]
        for(j in 1:ncol(phenoDat)){
            phenoDat[,j] = as.factor(as.character(phenoDat[,j]))
        }
        splitObject@phenoData = phenoDat
        splitObject@rt$raw = xcmsObject@rt$corrected[split.samples]
        splitObject@rt$corrected = xcmsObject@rt$corrected[split.samples]
        splitObject@rt$Original_raw = xcmsObject@rt$raw[split.samples]
        splitObject@filepaths = xcmsObject@filepaths[split.samples]
        
        xcmsObject.split[[k]] = splitObject
        #mutual.RT.correction[[k]] = list(raw = xcmsObject@rt$raw[split.samples], corrected = xcmsObject@rt$corrected[split.samples])
    }
    
    names(xcmsObject.split) = paste("xcms_",BatchNames,sep ="")
    #names(mutual.RT.correction) = paste("xcms_",BatchNames,sep ="")
    
    return(xcmsObject.split)
    
}

