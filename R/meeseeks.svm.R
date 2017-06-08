Meeseeks.SVM = function(FeatureMatrix, GroupLabels, SampleLabels = NULL, nFolds = 10, nSims = 20, plot.out = TRUE, nCPU = -1, plotcol = NULL, svm.kernel = "linear"){
    
    
    
    #FeatureMatrix = BreastCancer[,2:10]
    
    #GroupLabels = BreastCancer$Class
    # SampleLabels = BreastCancer$Id
    
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
        plotcol = "SVM"
    }
    
    if(!"matrix" %in% FeatureMatrix){
        warning("'FeatureMatrix' was not a matrix. Conversion applied.")
        FeatureMatrix = data.matrix(FeatureMatrix)
    }
    
    removed.entries = FALSE
    if(anyNA(FeatureMatrix)){
        warning("No missing values allowed. The non complete cases will be ignored.", immediate. = FALSE)
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
    #pb <- txtProgressBar(max=nSamp, style=3)
    #progress <- function(n) setTxtProgressBar(pb, n)
    #opts <- list(progress=progress)
    
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
            
            for (i in 1:nFolds){
                
                validData <- FeatureMatrix[Group.fold[[i]],]
                validLabels <- GroupLabels[Group.fold[[i]]]
                trainData <- FeatureMatrix[-Group.fold[[i]],]
                trainLabels <- GroupLabels[-Group.fold[[i]]]
                
                meeseeks.model <- e1071::svm(x = trainData, y= trainLabels, type = "C-classification", kernel = svm.kernel,probability = TRUE)
                
                predlabel.probs <- stats::predict(object = meeseeks.model, newdata = validData, decision.values = TRUE, probability = TRUE)
                
                #predlabel <- stats::predict(object = meeseeks.model, newdata = validData, type = "class")
                predcv.list[[i]] <- as.data.frame(cbind(SampleLabels[Group.fold[[i]]], attr(predlabel.probs, "probabilities")[,2,drop=FALSE])) 
                names(predcv.list[[i]])<-c("ID","decision.values")
                #predcv.votes.list[[i]] <- as.data.frame(cbind(rownames(validData), as.character(predlabel)))
                
                
                
            }
            
            
            
            # evaluate 1 run performance
            
            predcv.df <- data.table::rbindlist(predcv.list)
            ord <-order(as.numeric(unlist(Group.fold)))
            predcv.df<-predcv.df[ord,]
            
            
            
            predcv.df$decision.values<- as.numeric(as.character(predcv.df$decision.values))
            
            rocpred <- ROCR::prediction(predictions = predcv.df$decision.values, labels=GroupLabels)
            
            perf.roc <- ROCR::performance(rocpred, measure = "tpr", x.measure = "fpr")
            #perf.pr <- ROCR::performance(rocpred, measure = "prec", x.measure = "rec")
            
            perf.roc.xy=data.frame(perf.roc@x.values,perf.roc@y.values)
            colnames(perf.roc.xy)=c("minspecificity","sensitivity")
            random_diag_line=data.frame(c(0,1),c(0,1))
            colnames(random_diag_line)=c("x","y")
            
            #perf.pr.xy=data.frame(perf.pr@x.values,perf.pr@y.values)
            #colnames(perf.pr.xy)=c("precision","recall")
            
            #ggplot(perf.roc.xy, aes(x=minspecificity, y=sensitivity)) +geom_line(colour="red") +geom_line(data=random_diag_line, aes(x=x,y=y), colour="black")
            #ggplot(perf.pr.xy, aes(x=precision, y=recall)) +geom_line(colour="red") 
            
            aucvalue <- ROCR::performance(rocpred , measure="auc")
            
            perf.roc.xy$AUC = aucvalue@y.values[[1]]
            
            results[[iPar]] =perf.roc.xy
            
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
    
    perflist = list()
    AUCs = rep(NA, nSims)
    
    cter = 1
    for(main.loop in 1:length(performanceList)){
        for(small.loop in 1:length(performanceList[[main.loop]])){
            ROCdata = performanceList[[main.loop]][[small.loop]][,1:2]
            ROCdata$sim = cter
            perflist[[cter]] = ROCdata
            AUCs[cter] = performanceList[[main.loop]][[small.loop]][1,3]
            cter = cter + 1
        }
    }
    
    Performance = rbindlist(perflist)
    
    Nplotpoints = 100
    
    RC = data.frame(ROCx = seq(0, 1, length.out = Nplotpoints),
                    ROCy = rep(NA,Nplotpoints),
                    lower = rep(NA,Nplotpoints),
                    upper = rep(NA,Nplotpoints))
    
    
    
    for( l in 1: nrow(RC)){
        
        if(l!=1){
            ROC.data = Performance[Performance$minspecificity > RC$ROCx[l-1] & Performance$minspecificity <= RC$ROCx[l],  ]
        }else{
            ROC.data = Performance[Performance$minspecificity <= RC$ROCx[l],  ]
        }
        RC[l,2:4] = quantile(ROC.data$sensitivity, c(1/2,0.025,0.975))
        
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
    
    
    
    
    pp <- ggplot() + 
        geom_line(data = RC, aes(ROCx, ROCy,colour=plotcol)) + 
        geom_ribbon(data = RC, aes(x=ROCx, ymin=lower, ymax=upper, fill = plotcol), alpha=0.2) +
        geom_line(data=random_diag_line, aes(x=x,y=y),colour="black") +
        guides(colour=guide_legend(title="Method"), fill = guide_legend(title="95% interval"))+
        xlab("False Positive Rate (1-specificity)") +
        ylab("True Positive Rate (Sensitivity)") +
        ggtitle("ROC curve for SVM classification") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) 
    
    if(plot.out){
        print(pp)
    } 
    
    return(RC)
    
}