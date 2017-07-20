#' Restructure BatchProcessing and calculate combined results
#'
#'
#' @param PLSDAresults.list list of PLSDA variable importances
#' @param BatchProcessingStructure A BatchProcessingStructure
#' @param original.column.indices.list list with the original column indices
#' @param GroupData.list list with the original GroupData
#' @param FoldChanges.list list with fold change vectors
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
BatchResultsMerger = function(PLSDAresults.list, BatchProcessingStructure, original.column.indices.list, GroupData.list, FoldChanges.list){
    
    list.initiation = vector("list", length = nrow(BatchProcessingStructure))
    
    PLSDAresults.restructured             = list.initiation
    original.column.indices.restructured  = list.initiation
    GroupData.restructured                = list.initiation
    GroupData.restructured.relevantGroups = list.initiation
    FoldChanges.restructured              = list.initiation
    cter = 1
    for(l in 1:nrow(BatchProcessingStructure)){
        sublist.initiation = vector("list", length = BatchProcessingStructure$nBatches[l])
        
        PLSDAresults.restructured[[l]]             = sublist.initiation
        original.column.indices.restructured[[l]]  = sublist.initiation
        GroupData.restructured[[l]]                = sublist.initiation
        GroupData.restructured.relevantGroups[[l]] = sublist.initiation
        FoldChanges.restructured[[l]]              = sublist.initiation
        
        for(sublist in 1:BatchProcessingStructure$nBatches[l]){
            PLSDAresults.restructured[[l]][[sublist]] = PLSDAresults.list[[cter]]
            original.column.indices.restructured[[l]][[sublist]] = original.column.indices.list[[cter]]
            GroupData.restructured[[l]][[sublist]] = GroupData.list[[cter]]
            GroupData.restructured.relevantGroups[[l]][[sublist]] = GroupData.list[[cter]][original.column.indices.list[[cter]], ]
            FoldChanges.restructured[[l]][[sublist]] = FoldChanges.list[[cter]]
            cter = cter + 1
        }
        
    }
    
    PLSDA.interbatchResults = list.initiation
    for(l in 1:nrow(BatchProcessingStructure)){
        # improve on this for Nr Batch > 2 charlie
        
        Batch1.Groupdata  = as.data.frame(GroupData.restructured.relevantGroups[[l]][[1]][,1:6])
        groupMatches = vector("list", length = BatchProcessingStructure$nBatche - 1)
        groupMatches.notmatched = rep(NA,BatchProcessingStructure$nBatche - 1)
        for(sublist in 2:BatchProcessingStructure$nBatches[l]){
            if(BatchProcessingStructure$nBatches[l]>2){
                stop("this code only works for 2 batches so far")
            }
            BatchN.Groupdata  = as.data.frame(GroupData.restructured.relevantGroups[[l]][[sublist]][,1:6])    
            groupMatches[[sublist-1]] = MetaboMeeseeks::Std.matcher(Batch1.Groupdata$mzmed, Batch1.Groupdata$rtmed, as.matrix(BatchN.Groupdata[,c("mzmed","rtmed")]), max.ppm = 10,  max.RT.narrow = 3, max.RT.wide = 8, accept.ppm.diff = 5)
            groupMatches.notmatched[sublist-1] = nrow(BatchN.Groupdata) - length(unique(groupMatches[[sublist-1]]$Matched[!is.na(groupMatches[[sublist-1]]$Matched)]))
        }
        interbatchResults = data.frame(matrix(NA,nrow = nrow(Batch1.Groupdata) + sum(unlist(groupMatches.notmatched)), ncol = 6))
        colnames(interbatchResults) = c("Batch1.index", "Batch2.index", "Batch1.plsdaScore", "Batch2.plsdaScore", "Batch1.foldchange", "Batch2.foldchange")
        interbatchResults$Batch1.index[1:length(original.column.indices.restructured[[l]][[1]])] = original.column.indices.restructured[[l]][[1]]
        interbatchResults$Batch2.index[1:length(original.column.indices.restructured[[l]][[1]])] = original.column.indices.restructured[[l]][[2]][groupMatches[[1]]$Matched]
        interbatchResults$Batch2.index[(length(original.column.indices.restructured[[l]][[1]])+1):length(interbatchResults$Batch2.index)] = original.column.indices.restructured[[l]][[2]][ which(!seq(1:length(original.column.indices.restructured[[l]][[2]])) %in% unique(groupMatches[[sublist-1]]$Matched[!is.na(groupMatches[[sublist-1]]$Matched)]))  ]  
        
        interbatchResults$Batch1.plsdaScore[1:length(original.column.indices.restructured[[l]][[1]])] = PLSDAresults.restructured[[l]][[1]]
        interbatchResults$Batch2.plsdaScore[1:length(original.column.indices.restructured[[l]][[1]])] = PLSDAresults.restructured[[l]][[2]][groupMatches[[1]]$Matched]
        interbatchResults$Batch2.plsdaScore[(length(original.column.indices.restructured[[l]][[1]])+1):length(interbatchResults$Batch2.index)] = PLSDAresults.restructured[[l]][[2]][ which(!seq(1:length(original.column.indices.restructured[[l]][[2]])) %in% unique(groupMatches[[sublist-1]]$Matched[!is.na(groupMatches[[sublist-1]]$Matched)]))  ]
        
        interbatchResults$Batch1.foldchange[1:length(original.column.indices.restructured[[l]][[1]])] = as.numeric(FoldChanges.list[[l]][[1]])
        interbatchResults$Batch2.foldchange[1:length(original.column.indices.restructured[[l]][[1]])] = as.numeric(FoldChanges.list[[l]][[2]][groupMatches[[1]]$Matched])
        interbatchResults$Batch2.foldchange[(length(original.column.indices.restructured[[l]][[1]])+1):length(interbatchResults$Batch2.index)] = as.numeric(FoldChanges.list[[l]][[2]][ which(!seq(1:length(original.column.indices.restructured[[l]][[2]])) %in% unique(groupMatches[[sublist-1]]$Matched[!is.na(groupMatches[[sublist-1]]$Matched)]))  ])
        
        PLSDA.interbatchResults[[l]] = interbatchResults
    }
    
    return(PLSDA.interbatchResults)
    
}  


