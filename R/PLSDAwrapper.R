#' Function to split up xcms object by catch to perform a within batch analysis and merge the results at the end
#'
#' binnen batch analyse 
#' gebruik obiwarp om te weten hoe je de features van batch 1 en 2 aan elkaar moet plakken
#' doe dan analyse loess, PLSDA van batch 1 en 2 afzonderlijk. Giet alles in data frame (significant, fold change)
# 
#'
#' @param plsdaMatrix The matrix to be used inthe PLSDA
#' @param classLabels The class labels. The length must correpsond to the amount of rows in the plsdaMatrix
#' @param plotMarks (optional) A character vector with plot labels (same length as classLabels).
#' 
#'   
#' @return varimportance.PLSDA The vector with variable importances
#' 
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' 
#' 
#' @export
PLSDAwrapper <- function(plsdaMatrix, classLabels, plotMarks = NULL, plotOutputs = FALSE){
    
    if(lenth(classLabels) != nrow(plsdaMatrix)){
        stop("the number of classlabels is not equal to the amount of rows in the plsdaMatrix.")
    }
    
    classLabels = classlabels[[part]][ classlabels[[part]] %in% unlist(subclasses)]
    
    if(is.null(plotMarks)){
        plotMarks = classLabels
    }
    
    colnames(plsdaMatrix) = paste("feat",as.character(seq(1,ncol(plsdaMatrix))), sep = "" )
    colnames(plsdaMatrix) = as.character(seq(1,ncol(plsdaMatrix)))
    
    #CAM.plsda.perf <- mixOmics::plsda(plsdaMatrix, as.factor(classLabels), ncomp = 4)
    #perf.plsda <- mixOmics::perf(CAM.plsda.perf, validation = 'Mfold', folds = 5,
    #                             progressBar = FALSE, nrepeat = 10)
    
    #plot(perf.plsda, overlay = 'measure', sd=TRUE)
    
    CAM.plsda <- mixOmics::plsda(plsdaMatrix, as.factor(classLabels), ncomp = 4)
    
    if(plotOutputs){
    mixOmics::plotIndiv(CAM.plsda , comp = c(1,2),
                        group = as.factor(classLabels), ind.names = plotMarks, 
                        ellipse = TRUE, legend = TRUE, title = ' PLSDA comp 1 - 2')
    
    mixOmics::plotLoadings(CAM.plsda, contrib = "max", comp = 1)
    }
    
    varimportance.PLSDA = abs(CAM.plsda$loadings$X[,1])
  
    return(varimportance.PLSDA)
    
}