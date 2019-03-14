#' Easy Cross-Validated Random Forest Binary Classification
#'
#' This function quickly performs a cross-validated Random Forest classification on a data matrix.
#'
#'  
#' @param FeatureMatrix The matrix of Features (obtained by using the xcms::groupval function). Matrix has to have columns for features and rows for samples. 
#' @param GroupLabels The group labels. If not a factor a conversion will be applied. 
#' @param SampleLabels (optional) unique sample identifier.
#' @param nFolds Number of cross validation folds.
#' @param nSims Number of simulations (every simulation has different folds)
#' @param plot.out Whether to print the ROC curve (default is TRUE). 
#' @param plot.type Type of plot output. "ROC" for receiver operator characteristic (default) or "PR" for precision-recall.
#' @param nCPU The number of cores to use (default is the maximum amount available minus 2)
#' @param plotcol (optional) colour to use for the plot
#' @param plottitle.extra Optional extra character string to be added to every plot title.
#' @param Nplotpoints The amount of points used to construct the plot.
#' @param ... Extra paremeters to be passed along to \code{\link[randomForest]{randomForest}}
#'  
#' @return 
#' A ROC plot (if plot.out = TRUE) and a list with 2 elements: 1) a data frame with the ROC plot data and 2) a matrix with the variable importance for each cross validated simulation (nFolds * nSims times). 
#' 
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#'
#' @importFrom foreach %dopar%
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach
#' @importFrom caret createFolds
#' @importFrom randomForest randomForest
#' @importFrom stats predict quantile complete.cases
#' @importFrom data.table rbindlist
#' @importFrom ROCR prediction performance
#' 
#' @import parallel
#' @import ggplot2
#'  
#' @export
Meeseeks.RF = function(FeatureMatrix, GroupLabels, SampleLabels = NULL, nFolds = 10, nSims = 20, plot.out = TRUE, plot.type = "ROC", nCPU = -1, plotcol = NULL, plottitle.extra = NULL, Nplotpoints = 501, ...){

    if(length(unique(GroupLabels)) != 2){
        stop("This function is only available in case of two groups in the data.")
    }
    
    if(nrow(FeatureMatrix) != length(GroupLabels)){
        stop("The amount of rows in FeatureMatrix is not identical to the length of GroupLabels. Note the samples must be in rows and the features in columns.")
    }
    
    if (nCPU == -1) {
        nCPU <- max(c(1,(parallel::detectCores(all.tests = FALSE, logical = TRUE) - 2)))
        nCPU <- min(nSims, nCPU)
    }
    
    if(!is.factor(GroupLabels)){
        GroupLabels = as.factor(as.character(GroupLabels))
    }
    
    if(min(as.numeric(table(as.numeric(GroupLabels)))) < nFolds){
        stop("The number of folds is greater than the amount of unique entries for one group.")
    }
    
    if(is.null(SampleLabels)){
        SampleLabels = paste("Sample",as.character(seq(1,nrow(FeatureMatrix))), sep="")
    } else if(length(unique(SampleLabels)) != length(SampleLabels)){
        warning("The SampleLabels argument cannot contain identical entries. Every sample identifier has to be unique.")
        SampleLabels = paste("Sample",as.character(seq(1,nrow(FeatureMatrix))))
    }
    
    if(is.null(plotcol)){
        plotcol = "Random Forest"
    }
    
    if(!"matrix" %in% class(FeatureMatrix)){
        warning("'FeatureMatrix' was not a matrix. Conversion applied.")
        FeatureMatrix = data.matrix(FeatureMatrix)
    }
    
    removed.entries = FALSE
    if(anyNA(FeatureMatrix)){
        warning("No missing values are allowed for Random Forest classification. These non complete cases will be ignored.", immediate. = FALSE)
        SampleLabels = SampleLabels[complete.cases(FeatureMatrix)]
        GroupLabels = GroupLabels[complete.cases(FeatureMatrix)]
        FeatureMatrix = FeatureMatrix[complete.cases(FeatureMatrix), ]
        removed.entries = TRUE
    }
    
    if(removed.entries){
        Sys.sleep(1)
    }
    
    Groups = unique(GroupLabels)
    
    # initialize parallel computations
    cl <- parallel::makeCluster(nCPU)
    doSNOW::registerDoSNOW(cl)
    performanceList <- list()
    parCounter <- NULL
    
    # divide the number of runs for each core
    run.division <- split(seq(1,nSims), seq(1,nSims) %% nCPU)
    
    performanceList <- foreach::foreach(parCounter = 1:nCPU,  .inorder = FALSE, .packages = c("caret", "e1071","stats","randomForest","ROCR")) %dopar%
    {
        results = list()
        for(iPar in 1:length(run.division[[parCounter]])){
            
            Group.fold <- caret::createFolds(GroupLabels, k = nFolds, returnTrain = FALSE)
            
            # the model 
            predcv.list <- vector("list",nFolds)
            #predcv.votes.list = vector("list",nFolds)
            meeseeks.varImportance <- matrix(NA,ncol = ncol(FeatureMatrix), nrow = nFolds)
            
            
            for (i in 1:nFolds){
                
                validData <- FeatureMatrix[Group.fold[[i]],]
                validLabels <- GroupLabels[Group.fold[[i]]]
                trainData <- FeatureMatrix[-Group.fold[[i]],]
                trainLabels <- GroupLabels[-Group.fold[[i]]]
                
                meeseeks.model <- randomForest::randomForest(x = trainData, y= trainLabels, importance = TRUE, ...)
                
                predlabel.probs <- stats::predict(object = meeseeks.model, newdata = validData, type="prob")
                
                predcv.list[[i]] <- as.data.frame(cbind(SampleLabels[Group.fold[[i]]], predlabel.probs[,2,drop=FALSE])) 
                names(predcv.list[[i]])<-c("ID","decision.values")

                meeseeks.varImportance[i,] = as.numeric(randomForest::importance(meeseeks.model, type = 1) )
                
            }
        
            # evaluate 1 run performance
            
            predcv.df <- data.table::rbindlist(predcv.list)
            ord <-order(as.numeric(unlist(Group.fold)))
            predcv.df<-predcv.df[ord,]
          
            predcv.df$decision.values<- as.numeric(as.character(predcv.df$decision.values))
            
            rocpred <- ROCR::prediction(predictions = predcv.df$decision.values, labels=GroupLabels)
            perf.roc <- ROCR::performance(rocpred, measure = "tpr", x.measure = "fpr")
            perf.pr <- ROCR::performance(rocpred, measure = "prec", x.measure = "rec")
            
            perf.roc.xy=data.frame(perf.roc@x.values,perf.roc@y.values)
            perf.pr.xy=data.frame(perf.pr@x.values,perf.pr@y.values)
            
            colnames(perf.roc.xy)=c("minspecificity","sensitivity")
            colnames(perf.pr.xy)=c("recall", "precision")
            
            if(is.nan(perf.pr.xy$precision[1])){
                perf.pr.xy$precision[1] = 1
            }
            
            random_diag_line=data.frame(c(0,1),c(0,1))
            colnames(random_diag_line)=c("x","y")

            aucvalue <- ROCR::performance(rocpred , measure="auc")
            
            perf.roc.xy$AUC = aucvalue@y.values[[1]]
            
            
            ROC.and.importance = list(ROC = perf.roc.xy, PR = perf.pr.xy, varImp = meeseeks.varImportance)
            
            results[[iPar]] =ROC.and.importance
            
        }    
        
        return(results)
    }
    #close(pb)
    parallel::stopCluster(cl)
    # output = list of lists
    
    
    ##########################################################################################
    ##########################################################################################
    ###########################                                ###############################
    ###########################             Plotting           ###############################
    ###########################                                ###############################
    ##########################################################################################
    ##########################################################################################
    
    perflistROC = list()
    perflistPR = list()
    AUCs = rep(NA, nSims)
    varImportance = matrix(NA, ncol = ncol(FeatureMatrix), nrow = nFolds*nSims)
    
    cter = 1
    if(purrr::vec_depth(performanceList) == 5){
        for(main.loop in 1:length(performanceList)){
            for(small.loop in 1:length(performanceList[[main.loop]])){
                ROCdata = performanceList[[main.loop]][[small.loop]]$ROC[,1:2]
                ROCdata_expanded = data.frame(minspecificity = seq(0, 1, length.out = Nplotpoints),
                                              sensitivity = rep(NA,Nplotpoints),
                                              sim = cter)
                
                for( l in 1: nrow(ROCdata_expanded)){
                    if(l!=1){
                        ROCdata.expanded = ROCdata[ROCdata$minspecificity > ROCdata_expanded$minspecificity[l-1] & ROCdata$minspecificity <= ROCdata_expanded$minspecificity[l], ]   
                    }else{
                        ROCdata.expanded = ROCdata[ ROCdata$minspecificity <= ROCdata_expanded$minspecificity[l], ]
                    }
                    if(nrow(ROCdata.expanded) != 0){
                        ROCdata_expanded$sensitivity[l] = mean(ROCdata.expanded$sensitivity)
                    } else{
                        ROCdata_expanded$sensitivity[l] = max(ROCdata$sensitivity[ ROCdata$minspecificity <= ROCdata_expanded$minspecificity[l]])
                    }
                }
                
                perflistROC[[cter]] = ROCdata_expanded
                
                PRdata = performanceList[[main.loop]][[small.loop]]$PR[,1:2] # is this true charlie? 4,5? 
                PRdata$sim = cter
                perflistPR[[cter]] = PRdata
                
                AUCs[cter] = performanceList[[main.loop]][[small.loop]]$ROC[1,3]
                varImportance[((cter-1)*nFolds + 1):(cter*nFolds),] = performanceList[[main.loop]][[small.loop]]$varImp
                cter = cter + 1
            }
        }
    } else if(purrr::vec_depth(performanceList) == 4){
        for(main.loop in 1:length(performanceList)){
            
                ROCdata = performanceList[[main.loop]]$ROC[,1:2]
                ROCdata_expanded = data.frame(minspecificity = seq(0, 1, length.out = Nplotpoints),
                                              sensitivity = rep(NA,Nplotpoints),
                                              sim = cter)
                
                for( l in 1: nrow(ROCdata_expanded)){
                    if(l!=1){
                        ROCdata.expanded = ROCdata[ROCdata$minspecificity > ROCdata_expanded$minspecificity[l-1] & ROCdata$minspecificity <= ROCdata_expanded$minspecificity[l], ]   
                    }else{
                        ROCdata.expanded = ROCdata[ ROCdata$minspecificity <= ROCdata_expanded$minspecificity[l], ]
                    }
                    if(nrow(ROCdata.expanded) != 0){
                        ROCdata_expanded$sensitivity[l] = mean(ROCdata.expanded$sensitivity)
                    } else{
                        ROCdata_expanded$sensitivity[l] =  max(ROCdata$sensitivity[ ROCdata$minspecificity <= ROCdata_expanded$minspecificity[l]])
                    }
                }
                
                perflistROC[[cter]] = ROCdata_expanded
                
                PRdata = performanceList[[main.loop]]$PR[,1:2] # is this true charlie? 4,5? 
                PRdata$sim = cter
                perflistPR[[cter]] = PRdata
                
                AUCs[cter] = performanceList[[main.loop]]$ROC[1,3]
                varImportance[((cter-1)*nFolds + 1):(cter*nFolds),] = performanceList[[main.loop]]$varImp
                cter = cter + 1
            
        }
    }else{
        stop("Parallel loop indexing failed.")
    }
    
    
    PerformanceROC = data.table::rbindlist(perflistROC)
    PerformancePR = data.table::rbindlist(perflistPR)

    
    RC = data.frame(ROCx = seq(0, 1, length.out = Nplotpoints),
                    ROCy = rep(NA,Nplotpoints),
                    lower = rep(NA,Nplotpoints),
                    upper = rep(NA,Nplotpoints))
    
    PR = data.frame(PRx = seq(0, 1, length.out = Nplotpoints),
                    PRy = rep(NA,Nplotpoints),
                    lower = rep(NA,Nplotpoints),
                    upper = rep(NA,Nplotpoints))
    
    for( l in 1: nrow(RC)){
        if(l!=1){
            ROC.data = PerformanceROC[PerformanceROC$minspecificity > RC$ROCx[l-1] & PerformanceROC$minspecificity <= RC$ROCx[l],  ]
        }else{
            ROC.data = PerformanceROC[PerformanceROC$minspecificity <= RC$ROCx[l],  ]
        }
        if(nrow(ROC.data) != 0){
            RC[l,2:4] = quantile(ROC.data$sensitivity, c(1/2,0.025,0.975))
        } else{
            RC[l,2:4] = RC[(l-1),2:4]
        }
    }
    RC = rbind(c(0,0,0,RC$upper[1]),RC)
    
    if(anyNA(RC)){
        incomplete.cases = which(!complete.cases(RC))
        for(cc in 1:length(incomplete.cases)){
            for(incomplete.col in which(is.na(RC[incomplete.cases[cc],]))){
                RC[incomplete.cases[cc],incomplete.col] = RC[incomplete.cases[cc]-1,incomplete.col]
            }
        }
    }
    
    for( l in 1: nrow(PR)){
        if(l!=1){
            PR.data = PerformancePR[PerformancePR$recall > PR$PRx[l-1] & PerformancePR$recall <= PR$PRx[l],  ]
        }else{
            PR.data = PerformancePR[PerformancePR$recall <= PR$PRx[l],  ]
        }
        PR[l,2:4] = quantile(PR.data$precision, c(1/2,0.025,0.975), na.rm = T)
    }
    
    if(anyNA(PR)){
        incomplete.cases = which(!complete.cases(PR))
        for(cc in 1:length(incomplete.cases)){
            for(incomplete.col in which(is.na(PR[incomplete.cases[cc],]))){
                PR[incomplete.cases[cc],incomplete.col] = PR[incomplete.cases[cc]-1,incomplete.col]
            }
        }
    }
    
    random_diag_line=data.frame(c(0,1),c(0,1))
    colnames(random_diag_line)=c("x","y")
    
    ppROC <- ggplot() + 
        geom_line(data = RC, aes_string(x = 'ROCx', y = 'ROCy', colour = 'plotcol')) + 
        geom_ribbon(data = RC, aes_string(x = 'ROCx', ymin = 'lower', ymax = 'upper', fill = 'plotcol'), alpha=0.2) +
        geom_line(data = random_diag_line, aes_string(x = 'x', y = 'y'),colour="black") +
        guides(colour=guide_legend(title="Method"), fill = guide_legend(title="95% interval"))+
        xlab("False Positive Rate (1-specificity)") +
        ylab("True Positive Rate (Sensitivity)") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_x_continuous(limits = c(0,1), expand = c(0.005,0.005)) +
        scale_y_continuous(limits = c(0,1), expand = c(0.005,0.005))
    if(is.null(plottitle.extra)){
        ppROC <- ppROC + ggtitle(paste("Random Forest ROC. Mean AUC = ", mean(AUCs),sep = "")) 
    }else{
        ppROC <- ppROC + ggtitle(paste("Random Forest ROC. Mean AUC = ", mean(AUCs),", ",plottitle.extra,sep = ""))
    }
    
    ppPR <- ggplot() + 
        geom_line(data = PR, aes_string(x = 'PRx', y = 'PRy', colour = 'plotcol')) + 
        geom_ribbon(data = PR, aes_string(x = 'PRx', ymin = 'lower', ymax = 'upper', fill = 'plotcol'), alpha=0.2) +
        guides(colour = guide_legend(title="Method"), fill = guide_legend(title="95% interval"))+
        xlab("Recall") +
        ylab("Precision") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_x_continuous(limits = c(0,1), expand = c(0.005,0.005)) +
        scale_y_continuous(limits = c(0,1), expand = c(0.005,0.005))
    if(is.null(plottitle.extra)){
        ppPR <- ppPR + ggtitle("Random Forest PR") 
    }else{
        ppPR <- ppPR + ggtitle(paste("Random Forest PR, ",plottitle.extra,sep = ""))
    }
    
    
    if(plot.out){
        if( plot.type != "PR"){
            print(ppROC)
        } else{
            print(ppPR)
        }
        
    } 
    
    Results = list( ROCdata = RC, varImportance = varImportance, ROCplot = ppROC, PRplor = ppPR, AUCs = AUCs)
    return(Results)
}

