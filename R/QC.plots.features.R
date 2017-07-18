#' QC plots for feature data
#'
#' Mean number of peaks are plotted for the peak data
#'
#'  
#' @param FeatureMatrix the matrix of Features (obtained by using the xcms::groupval function). Matrix has to have columns for features and rows for samples. 
#' @param XCMSobject An xcmsSet object
#' @param className In case there are multiple class definitions in the XCMSobject@phenoData object, the one of interest can be specified here. If not the code will ask.
#' @param NA.numeric.limit If missing values are not indicated as NA but with a numeric value like 0, supply the numerical value under which values are considered missing.
#'  
#' @return 
#' several QC plots
#' 
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#'
#'
#' @import ggplot2
#'  
#' @export
QC.plots.features = function(FeatureMatrix, XCMSobject, className = NULL, NA.numeric.limit = NULL){
 
    if(ncol(XCMSobject@phenoData) > 1 & is.null(className)){
        print(head(XCMSobject@phenoData))
        class <- readline(paste("Which of the following classes describes the individual groups: ",paste(colnames(XCMSobject@phenoData), collapse=" or ")))
    } else if (ncol(XCMSobject@phenoData) > 1 & !is.null(className)){
        if(className %in% colnames(XCMSobject@phenoData)){
            class = className
        } else{
            warning("The 'className' variable does not match with any names in the XCMSobject@phenoData object.")
            print(head(XCMSobject@phenoData))
            class <- readline(paste("Which of the following classes describes the individual groups: ",paste(colnames(XCMSobject@phenoData), collapse=" or ")))
        }
    } else { # only 1 available class
        class = colnames(XCMSobject@phenoData)
    }
    
    classes = as.character(unique(XCMSobject@phenoData[class][,1]))
    
    max.class.size = max(table(stringr::str_count(as.character(XCMSobject@phenoData[class][,1]))))
    missing.vals.df = expand.grid(class = classes, Nmissing = seq(0, max.class.size, 1),count = 0)
    colnames(missing.vals.df)[3] = "count"
    RSD.matrix = matrix(NA,ncol=ncol(FeatureMatrix), nrow = length(classes))
    no.missing.matrix = matrix(0,ncol=ncol(FeatureMatrix), nrow = length(classes))
    for(cl in 1:length(classes)){
        sample.nrs = which(XCMSobject@phenoData[class] == classes[cl])
        missing.values = rep(0,ncol(FeatureMatrix))
        for(ft in 1:ncol(FeatureMatrix)){
            missing.values[ft] = sum( FeatureMatrix[sample.nrs, ft] < 0.1)
            # if at least 80 percent of samples are available, calculate rsd
            if( sum(FeatureMatrix[sample.nrs, ft] < 0.1)/length(sample.nrs) <= 0.2){
                featVals = FeatureMatrix[sample.nrs, ft]
                RSD.matrix[cl,ft] = 100*sd(featVals[featVals> 0])/mean(featVals[featVals> 0])
            }
            if(sum(FeatureMatrix[sample.nrs, ft] < 0.1) < 1 ){
                no.missing.matrix[cl,ft] = 1
            }
            
        }
        for(nmis in 0:max.class.size){
            missing.vals.df$count[missing.vals.df$class == classes[cl] &missing.vals.df$Nmissing == nmis] = sum(missing.values == nmis)
        }
        
    }
    
    missing.vals.df2 = missing.vals.df
    missing.vals.df2$class = as.character(missing.vals.df2$class)
    #missing.vals.df2 = rbind(missing.vals.df2, c(" common", 0, NA),c(" common", 1, NA),c(" common", 2, NA),c(" common", 3, NA),c(" common", 4, NA),c(" common", 5, NA),c(" common", 6, NA))
    missing.vals.df2$class = as.factor(missing.vals.df2$class)
    missing.vals.df2$Nmissing= as.numeric(missing.vals.df2$Nmissing)
    missing.vals.df2$count = as.numeric(missing.vals.df2$count)
    
    common.peaks.none.missing = sum(colSums( no.missing.matrix) == 5)
    common.sub.df = missing.vals.df2[missing.vals.df2$Nmissing==0,]
    common.sub.df$common = "common"
    common.sub.df$count[common.sub.df$class != " common"] = common.peaks.none.missing
    gg1 <- ggplot() +  
        geom_bar(data = missing.vals.df2, 
                 aes(x = Nmissing, y = count, fill = class),
                 position="dodge", 
                 stat="identity", 
                 width = 0.6,
                 show.legend = TRUE) +
        theme_bw() +
        geom_bar(data = common.sub.df, 
                 aes(x = Nmissing, y = count,  fill = class, colour = common), 
                 position="dodge", 
                 stat="identity", 
                 width = 0.6, 
                 colour = "black",
                 show.legend = FALSE) 
        
    print(gg1)
    
 
    RSD.plotdata = data.frame(matrix(NA,nrow=0,ncol=2))
    
    for(cl in 1:length(classes)){
        RSDs = RSD.matrix[cl,!is.na(RSD.matrix[cl,])]
        RSD.plotdata = rbind(RSD.plotdata, cbind(rep(as.character(classes[cl]),length(RSDs)),RSDs))
        
    }
    colnames(RSD.plotdata)=c("class","RSDs")
    RSD.plotdata$RSDs = as.numeric(as.character(RSD.plotdata$RSDs))
    gg2 <- ggplot(RSD.plotdata, aes(x=RSDs, colour=as.factor(class))) + 
        geom_density(adjust = 1) +
        theme_bw() +
        geom_hline(color = "#e6e6e6", yintercept = 0, size= 0.7)
    # axis ticks per 20, median RSD after legend
    
    print(gg2)
    
    Npeaks_pr_class = rep(NA,length(classes))
    for(cls in 1:length(classes)){
        Npeaks_pr_class[cls] = nrow(RSD.plotdata[RSD.plotdata$class == classes[cls],])
    }
    npeaks.df = data.frame(class = classes, Npeaks = Npeaks_pr_class)
    gg3 <- ggplot(RSD.plotdata, aes(y=RSDs, x =class)) + 
        geom_boxplot() +
        geom_text(data= npeaks.df, aes(x=classes, y=0, label=Npeaks),colour = "red") +
        theme_bw()
    
    print(gg3)
    
    mRSD <- data.frame(cbind(as.character(classes),rowMedians(RSD.matrix,na.rm=TRUE)))
    colnames(mRSD) = c("class","median RSD")
    print(mRSD)
    
}

