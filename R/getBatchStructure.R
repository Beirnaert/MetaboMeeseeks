#' Get the Batch processing Structure. Necessary for merging batches again after analysis
#'
#'
#' @param xcmsSets.Post.batchSplit list of xcmsSet objects after Batch splitting (BatchSplitR function)
#' 
#' @return BatchProcessingStructure A data frame with the list structure (quantifies the list nesting)
#' 
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' 
#' 
#' @export
#' 
#' @importFrom purrr map
getBatchStructure = function(xcmsSets.Post.batchSplit){
    
    BatchProcessingStructure = data.frame(MultiBatchData = 1:length(xcmsSets.Post.batchSplit), nBatches = NA)
    BatchProcessingStructure$nBatches = as.numeric(
        unlist(
            purrr::map(.x = xcmsSets.Post.batchSplit,
                       .f = ~ length(.x)
            )
        )
    )
    
    return(BatchProcessingStructure)
}