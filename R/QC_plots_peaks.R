#' QC plots for peak data
#'
#' Mean number of peaks are plotted for the peak data
#'
#'  
#' @param XCMSobject An xcmsSet object
#' @param include.postfill.plots Whether to include the peaks added after filling (only if available)
#' @param className In case there are multiple class definitions in the XCMSobject@phenoData object, the one of interest can be specified here. If not the code will ask.
#' @param plottitle (optional) The title to place above the plot(s). 
#'  
#' @return 
#' A QC plot for nr of peaks
#' 
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#'
#'
#' @import ggplot2
#' @importFrom utils head
#' @importFrom stats sd
#'  
#' @export

QC_plots_peaks <- function(XCMSobject, include_postfill_plots = FALSE, className = NULL, plottitle = NULL){
    
    if(ncol(XCMSobject@phenoData) > 1 & is.null(className)){
        print(head(XCMSobject@phenoData))
        class <- readline(paste("Which of the following classes describes the individual groups: ",paste(colnames(XCMSobject@phenoData), collapse=" or ")))
    } else if (ncol(XCMSobject@phenoData) > 1 & !is.null(className)){
        if(className %in% colnames(XCMSobject@phenoData)){
            class <- className
        } else{
            warning("The 'className' variable does not match with any names in the XCMSobject@phenoData object.")
            print(head(XCMSobject@phenoData))
            class <- readline(paste("Which of the following classes describes the individual groups: ",paste(colnames(XCMSobject@phenoData), collapse=" or ")))
        }
    } else { # only 1 available class
        class <- colnames(XCMSobject@phenoData)
    }
    
    classes <- as.character(unique(XCMSobject@phenoData[class][,1]))
    class_meanNpeaks_prefil <- rep(NA,length(classes))
    class_sdNpeaks_prefil <- rep(NA,length(classes))
    class_meanNpeaks_postfil <- rep(NA,length(classes))
    class_sdNpeaks_postfil <- rep(NA,length(classes))
    peaks_df <- data.frame(XCMSobject@peaks)
    for(cl in 1:length(classes)){
        sample_nrs <- which(XCMSobject@phenoData[class] == classes[cl])
        smpl_npeaks_prefil <- rep(NA,length(sample_nrs))
        smpl_npeaks_postfil <- rep(NA,length(sample_nrs))
        for(smpl in 1:length(sample_nrs)){
            smpl_npeaks_prefil[smpl] <- nrow(peaks_df[peaks_df$sample==sample_nrs[smpl] & !is.na(peaks_df$intb),])
            if(include_postfill_plots){
                smpl_npeaks_postfil[smpl] <- nrow(peaks_df[peaks_df$sample==sample_nrs[smpl] &  peaks_df$into>0,])
            }
        }
        class_meanNpeaks_prefil[cl] <- mean(smpl_npeaks_prefil)
        class_sdNpeaks_prefil[cl] <- sd(smpl_npeaks_prefil)
        if(include_postfill_plots){
            class_meanNpeaks_postfil[cl] <- mean(smpl_npeaks_postfil)
            class_sdNpeaks_postfil[cl] <- sd(smpl_npeaks_postfil)
        }
    }
    
    prefil_df <- data.frame(cbind(as.character(classes), class_meanNpeaks_prefil, class_sdNpeaks_prefil))
    colnames(prefil_df) <- c("class","mean","sd")
    prefil_df$mean <- as.numeric(as.character(prefil_df$mean))
    prefil_df$sd <- as.numeric(as.character(prefil_df$sd))
    g1 <- ggplot(prefil_df, aes(x = class, y = mean)) +  
          geom_bar(position="dodge", stat="identity", width = 0.6,color="black", fill = "grey") + 
          geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) + 
          theme_bw() +   
          scale_y_continuous(name="Mean number of peaks (+/- sd) (no peak filling)")  
    if(!is.null(plottitle)){
        g1 <- g1 + ggtitle(plottitle) +
            theme(plot.title = element_text(hjust = 0.5))
    }
    print(g1)
    
    if(include_postfill_plots){
        postfil_df <- data.frame(cbind(as.character(classes),class_meanNpeaks_postfil,class_sdNpeaks_postfil))
        colnames(postfil_df) <- c("class","mean","sd")
        postfil_df$mean <- as.numeric(as.character(postfil_df$mean))
        postfil_df$sd <- as.numeric(as.character(postfil_df$sd))
        g2 <- ggplot(prefil_df, aes(x = class, y = mean)) +  
              geom_bar(position="dodge", stat="identity", width = 0.6,color="black", fill = "grey") + 
              geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) + 
              theme_bw() +   
              scale_y_continuous(name="Mean number of peaks (+/- sd) (with peak filling)")  
        if(!is.null(plottitle)){
            g2 <- g2 + ggtitle(plottitle) +
                theme(plot.title = element_text(hjust = 0.5))
        }
        print(g2)
    }

    return(g1)
}

