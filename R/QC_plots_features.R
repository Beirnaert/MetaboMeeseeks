#' QC plots for feature data
#'
#' RSD's for each class are plotted.
#'
#'  
#' @param FeatureMatrix the matrix of Features (obtained by using the xcms::groupval function). Matrix has to have columns for features and rows for samples. 
#' @param classVector The vector of sample classes. 
#' @param COI (optional) In case only a specific class is of interest.
#' @param NA_numeric_limit If missing values are not indicated as NA but with a numeric value like 0, supply the numerical value under which values are considered missing.
#' @param plottitle (optional) Plot title 
#' @param BOI_vector (optional) In case the RSD's need only be calculated for a specific subset of batches, provide the character/factor vector witch batches here (same length as the number of rows in FeatureMatrix).
#' @param BOI (optional) The name of the batch of interest. BOI_vector has to be supplied as well
#' @param RDSplot_text_adjustheigt Height adjustment for the count number in the RSD boxplot visuzlization.
#'  
#' @return 
#' several QC plots
#' 
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#'
#'
#' @import ggplot2
#'  
#' @export
#' 
#' @importFrom utils head
#' @importFrom matrixStats rowMedians
#' 
QC_plots_features <- function(FeatureMatrix, classVector, COI = NULL, NA_numeric_limit = NULL, plottitle = NULL, BOI_vector = NULL, BOI = NULL, RDSplot_text_adjustheigt = 0){
 
    
    if(length(classVector) != nrow(FeatureMatrix)){
        stop("the length of classVector is not equal to the number of rows in FeatureMatrix.")
    }
    
    if(!is.null(COI)){
        if(COI %in% classVector){
            FeatureMatrix <- FeatureMatrix[classVector == COI,]
            if(!is.null(BOI_vector)){
                BOI_vector <- BOI_vector[classVector == COI]
            }
            classVector <- classVector[classVector == COI]
        }else{
            stop("COI is does not match with anything in classVector.")
        }
    }
    
    if(!is.null(NA_numeric_limit)){
        if(is.numeric(NA_numeric_limit)){
            FeatureMatrix[FeatureMatrix < NA_numeric_limit] <- NA
        }
    }
    
    batchspecific = ""
    if(!is.null(BOI_vector) & !is.null(BOI) ){
        if(length(BOI_vector) == nrow(FeatureMatrix)){
            batchspecific <- paste(",", BOI, sep = " ")
            FeatureMatrix <- FeatureMatrix[BOI_vector ==  BOI ,]
            classVector <- classVector[BOI_vector ==  BOI]
        }else{
            stop("length of supplied BOI_vector is not the same as the number of rows in FeatureMatrix.")
        }
    }
    
    classes <- as.character(unique(classVector))
    max_class_size <- max(summary(as.factor(classVector)))
    missing_vals_df <- expand.grid(class = classes, Nmissing = seq(0, max_class_size, 1),count = 0)
    colnames(missing_vals_df)[3] <- "count"
    RSD_matrix <- matrix(NA,ncol=ncol(FeatureMatrix), nrow = length(classes))
    no_missing_matrix <- matrix(0,ncol=ncol(FeatureMatrix), nrow = length(classes))
    
    for(cl in 1:length(classes)){
        sample.nrs <- which(classVector == classes[cl])
        missing.values <- rep(0,ncol(FeatureMatrix))
        for(ft in 1:ncol(FeatureMatrix)){
            missing.values[ft] <- sum(is.na(FeatureMatrix[sample.nrs, ft]))  #sum( FeatureMatrix[sample.nrs, ft] < 0.1)
            # if at least 80 percent of samples are available, calculate rsd
            if( missing.values[ft]/length(sample.nrs) <= 0.2){
                featVals <- as.numeric(na.omit(FeatureMatrix[sample.nrs, ft]))
                RSD_matrix[cl,ft] <- 100*sd(featVals)/mean(featVals)
            }
            if(missing.values[ft] < 1 ){
                no_missing_matrix[cl,ft] <- 1
            }
            
        }
        for(nmis in 0:max_class_size){
            missing_vals_df$count[missing_vals_df$class == classes[cl] &missing_vals_df$Nmissing == nmis] <- sum(missing.values == nmis)
        }
        
    }
    
    missing_vals_df2 <- missing_vals_df
    missing_vals_df2$class <- as.character(missing_vals_df2$class)
    #missing_vals_df2 = rbind(missing_vals_df2, c(" common", 0, NA),c(" common", 1, NA),c(" common", 2, NA),c(" common", 3, NA),c(" common", 4, NA),c(" common", 5, NA),c(" common", 6, NA))
    missing_vals_df2$class <- as.factor(missing_vals_df2$class)
    missing_vals_df2$Nmissing <- as.numeric(missing_vals_df2$Nmissing)
    missing_vals_df2$count <- as.numeric(missing_vals_df2$count)
    
    common_peaks_none_missing <- sum(colSums( no_missing_matrix) == 5)
    common_sub_df <- missing_vals_df2[missing_vals_df2$Nmissing==0,]
    common_sub_df$common <- "common"
    common_sub_df$count[common_sub_df$class != " common"] <- common_peaks_none_missing
    gg1 <- ggplot() +  
        geom_bar(data = missing_vals_df2, 
                 aes_string(x = 'Nmissing', y = 'count', fill = 'class'),
                 position="dodge", 
                 stat="identity", 
                 width = 0.6,
                 show.legend = TRUE) +
        theme_bw() +
        geom_bar(data = common_sub_df, 
                 aes_string(x = 'Nmissing', y = 'count',  fill = 'class', colour = 'common'), 
                 position="dodge", 
                 stat="identity", 
                 width = 0.6, 
                 colour = "black",
                 show.legend = FALSE) 
        
    if(!is.null(plottitle)){
        gg1 <- gg1 + ggtitle(paste(plottitle, batchspecific, sep = " ") ) +
            theme(plot.title = element_text(hjust = 0.5))
    }
    print(gg1)
    
 
    RSD_plotdata <- data.frame(matrix(NA,nrow=0,ncol=2))
    
    for(cl in 1:length(classes)){
        RSDs <- RSD_matrix[cl,!is.na(RSD_matrix[cl,])]
        RSD_plotdata <- rbind(RSD_plotdata, cbind(rep(as.character(classes[cl]),length(RSDs)),RSDs))
        
    }
    colnames(RSD_plotdata) <- c("class","RSDs")
    RSD_plotdata$RSDs <- as.numeric(as.character(RSD_plotdata$RSDs))
    gg2 <- ggplot(RSD_plotdata, aes_string(x = "RSDs", colour = "as.factor(class)")) + 
        geom_density(adjust = 1) +
        theme_bw() +
        geom_hline(color = "#e6e6e6", yintercept = 0, size= 0.7)
    # axis ticks per 20, median RSD after legend
    if(!is.null(plottitle)){
        gg2 <- gg2 + ggtitle( paste(plottitle, batchspecific, sep = " ")  ) +
            theme(plot.title = element_text(hjust = 0.5))
    }
    print(gg2)
    
    Npeaks_pr_class <- rep(NA,length(classes))
    for(cls in 1:length(classes)){
        Npeaks_pr_class[cls] <- nrow(RSD_plotdata[RSD_plotdata$class == classes[cls],])
    }
    npeaks_df <- data.frame(class = classes, Npeaks = Npeaks_pr_class)
    gg3 <- ggplot(RSD_plotdata, aes_string(y = "RSDs", x = "class")) + 
        geom_boxplot() +
        geom_text(data= npeaks_df, aes_string(x = "classes", y = 0+RDSplot_text_adjustheigt, label = "Npeaks"),colour = "red") +
        theme_bw()
    if(!is.null(plottitle)){
        gg3 <- gg3 + ggtitle( paste(plottitle, batchspecific, sep = " ")  ) +
            theme(plot.title = element_text(hjust = 0.5))
    }
    print(gg3)
    
    mRSD <- data.frame(cbind(as.character(classes),rowMedians(RSD_matrix,na.rm=TRUE)))
    colnames(mRSD) <- c("class","median RSD")
    print(mRSD)
    return(gg3)
}

