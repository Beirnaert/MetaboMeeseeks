#' QC plots for peak data
#'
#' Mean number of peaks are plotted for the peak data
#'
#'  
#' @param XCMSobject An xcmsSet object
#' @param include.postfill.plots Whether to include the peaks added after filling (only if available)
#' @param className In case there are multiple class definitions in the XCMSobject@phenoData object, the one of interest can be specified here. If not the code will ask.
#'  
#' @return 
#' a QC plot for nr of peaks
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

QC.plots.peaks = function(XCMSobject, include.postfill.plots = FALSE, className = NULL){
    
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
    class.meanNpeaks.prefil = rep(NA,length(classes))
    class.sdNpeaks.prefil = rep(NA,length(classes))
    class.meanNpeaks.postfil = rep(NA,length(classes))
    class.sdNpeaks.postfil = rep(NA,length(classes))
    peaks.df = data.frame(XCMSobject@peaks)
    for(cl in 1:length(classes)){
        sample.nrs = which(XCMSobject@phenoData[class] == classes[cl])
        smpl.npeaks.prefil = rep(NA,length(sample.nrs))
        smpl.npeaks.postfil = rep(NA,length(sample.nrs))
        for(smpl in 1:length(sample.nrs)){
            smpl.npeaks.prefil[smpl] = nrow(peaks.df[peaks.df$sample==sample.nrs[smpl] & !is.na(peaks.df$intb),])
            if(include.postfill.plots){
                smpl.npeaks.postfil[smpl] = nrow(peaks.df[peaks.df$sample==sample.nrs[smpl] &  peaks.df$into>0,])
            }
        }
        class.meanNpeaks.prefil[cl] = mean(smpl.npeaks.prefil)
        class.sdNpeaks.prefil[cl] = sd(smpl.npeaks.prefil)
        if(include.postfill.plots){
            class.meanNpeaks.postfil[cl] = mean(smpl.npeaks.postfil)
            class.sdNpeaks.postfil[cl] = sd(smpl.npeaks.postfil)
        }
    }
    
    prefil.df = data.frame(cbind(as.character(classes),class.meanNpeaks.prefil,class.sdNpeaks.prefil))
    colnames(prefil.df) = c("class","mean","sd")
    prefil.df$mean = as.numeric(as.character(prefil.df$mean))
    prefil.df$sd = as.numeric(as.character(prefil.df$sd))
    g1 <- ggplot(prefil.df, aes(x = class, y = mean)) +  
          geom_bar(position="dodge", stat="identity", width = 0.6,color="black", fill = "grey") + 
          geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) + 
          theme_bw() +   
          scale_y_continuous(name="Mean number of peaks (± sd) (no peak filling)")  
    print(g1)
    
    if(include.postfill.plots){
        postfil.df = data.frame(cbind(as.character(classes),class.meanNpeaks.postfil,class.sdNpeaks.postfil))
        colnames(postfil.df) = c("class","mean","sd")
        postfil.df$mean = as.numeric(as.character(postfil.df$mean))
        postfil.df$sd = as.numeric(as.character(postfil.df$sd))
        g2 <- ggplot(prefil.df, aes(x = class, y = mean)) +  
              geom_bar(position="dodge", stat="identity", width = 0.6,color="black", fill = "grey") + 
              geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) + 
              theme_bw() +   
              scale_y_continuous(name="Mean number of peaks (± sd) (with peak filling)")  
        print(g2)
    }

}

