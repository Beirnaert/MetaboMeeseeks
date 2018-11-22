#' Get EIC's and plot in ggplot
#'
#' 
#'
#' @param XCMSobject A single XCMS object
#' @param groups (optional) The group indices for which to obtain an EIC. Default is all groups.
#' @param sample_names (optional) The names of samples (from XCMSobject@phenoData) to include in the EIC plots.
#' @param sample_indices (optional) The indices of samples (from XCMSobject@phenoData) to include in the EIC plots.
#' @param phenotype_classname Variable name (from XCMSobject@phenoData) containing class information for the plot colours. 
#' @param plot_samples Whether to plot individual samples (default = FALSE: sample groups plotted in distinct colours)
#' @param ... Other variables to be passed along to \link[xcms]{getEIC}
#'
#'   
#' @return EICresults A list with the data for plotting, and the ggplot objects.
#' 
#'
#' @author Matthias Cuyckx, \email{matthias.cuykx@@uantwerpen.be}
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @importFrom xcms getEIC
#' @import ggplot2
#'  
#' 
#' @export
EICsBatchPlotter <- function(XCMSobject, groups = NULL, sample_names = NULL, sample_indices = NULL, phenotype_classname = NULL, plot_samples = FALSE, ...){
    
    ## checks for the arguments
    if(is.null(groups)){
        groups <- 1:length(XCMSobject@groupidx)
    } else{
        if(max(groups) > length(XCMSobject@groupidx)){
            stop("The 'groups' variable contains a group number larger than the amount of groups in the XCMSobject")
        }
    }
    
    if(is.null(sample_names) & is.null(sample_indices)){
        sample_indices <- seq(1,length(rownames(XCMSobject@phenoData)))
    } else if(is.null(sample_indices) & !is.null(sample_names)){
        sample_indices <- which(rownames(XCMSobject@phenoData) %in% sample_names)
        if(length(sample_names) != length(sample_indices)){
            warning("Certain sample_names could not be matched to the sample names in the XCMSobject (rownames of @phenoData)")
        }
    }
    
    if(max(sample_indices) > length(rownames(XCMSobject@phenoData))){
        stop("There are more 'samples' than the amount of samples in the XCMSobject(@phenoData)")
    }
    sample_indices <- as.integer(sample_indices)
    
    if(!plot_samples){
        if(!is.null(phenotype_classname) & !phenotype_classname %in% names(XCMSobject@phenoData)){
            phenotype_classname = names(XCMSobject@phenoData)[1]
            warning(paste("The 'phenotype_classname' variable is not in the phenoData of the XCMSobject. The variable '", phenotype_classname,"' is chosen", sep = ""))
        }else if(is.null(phenotype_classname)){
            phenotype_classname = names(XCMSobject@phenoData)[1]
            message(paste("The 'phenotype_classname' variable is chosen to be ",phenotype_classname, sep = ""))
        }
    }
    
    
    ## collecting EICs
    EICint <- xcms::getEIC(object = XCMSobject,
                           groupidx = groups,  
                           sampleidx = sample_indices, ...)
    
    
    ## grouping EICs across samples into dataframe for ggplot
    EICplots <- list()
    EICplotdata <- list()
    
    for (k in 1:length(EICint@eic[[1]])){
        eiclength <- NA
        for (l in 1:length(EICint@eic)){
            eiclength[l] <- dim(EICint@eic[[l]][[k]])[1]
        }
        
        EIC_df <- data.frame(rt = rep(NA,sum(eiclength)),
                             intensity = rep(NA,sum(eiclength)),
                             sample = rep(NA,sum(eiclength)),
                             class = rep(NA,sum(eiclength)))
        
        stop = 0
        for(m in 1:length(EICint@eic)){
            start <- stop +1
            stop <- sum(eiclength[1:m])
            
            if(!plot_samples){
                EIC_df$rt[start:stop] <- EICint@eic[[m]][[k]][,"rt"]
                EIC_df$intensity[start:stop] <- EICint@eic[[m]][[k]][,"intensity"]
                EIC_df$sample[start:stop] <- sample_indices[m]
                EIC_df$class[start:stop] <- as.character(XCMSobject@phenoData[sample_indices[m],phenotype_classname])
            } else{
                EIC_df$rt[start:stop] <- EICint@eic[[m]][[k]][,"rt"]
                EIC_df$intensity[start:stop] <- EICint@eic[[m]][[k]][,"intensity"]
                EIC_df$sample[start:stop] <- sample_indices[m]
            }
        }
        
        
        
        ## plotting
        if(!plot_samples){
            gg <- ggplot(EIC_df, aes_string(x = "rt", y = "intensity", color = "class", group = "sample")) + 
                geom_line() +
                ggtitle(paste("m/z",round(XCMSobject@groups[groups[k],"mzmed"],4),
                              "@",round(XCMSobject@groups[groups[k],"rtmed"],0), "s")) + 
                xlab("RT (s)") +
                ylab("raw intensities")
        }else{
            EIC_df$sample = as.factor(EIC_df$sample)
            gg <- ggplot(EIC_df, aes_string(x = "rt", y = "intensity", color = "sample", group = "sample")) + 
                geom_line() +
                ggtitle(paste("m/z",round(XCMSobject@groups[groups[k],"mzmed"],4),
                              "@",round(XCMSobject@groups[groups[k],"rtmed"],0), "s")) + 
                xlab("RT (s)") +
                ylab("raw intensities")
        }
        
        
        EICplots[[k]] <- gg
        EICplotdata[[k]] <- EIC_df
    }
    
    
    EICresults <- list(EICplots = EICplots, EICplotdata = EICplotdata)
    
    return(EICresults)
}
