#' obtain all datasets in a single list for easy processing with purrr
#'
#'
#' @param xcmsSets.Post.batchSplit list of xcmsSet objects after Batch splitting (BatchSplitR function)
#' @param BatchProcessingStructure The batch processing structure (obtained with getBatchStructure)
#' 
#' @return unStructuredBatchList a list with every individual batch in a list element
#' 
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' 
#' 
#' @export
unStructureBatches = function(xcmsSets.Post.batchSplit, BatchProcessingStructure){
    
    nMulitBatches = BatchProcessingStructure$MultiBatchData
    nBatches = sum(BatchProcessingStructure$nBatches)
    unStructuredBatchList = vector("list", length = nBatches )
    
    cter = 1
    for(k in seq_along(nMulitBatches)){
        for(j in 1:BatchProcessingStructure$nBatches[k]){
            unStructuredBatchList[[cter]] = xcmsSets.Post.batchSplit[[k]][[j]]
            cter = cter + 1
        }
    }
    
    return(unStructuredBatchList)
}