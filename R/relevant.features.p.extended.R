#' identify features (columns in the datamatrix) which are significantly associated with the outcome vector
#'
#' This function produces a p-value for every column in the datamatrix, corresponding to the null hypothesis that outcome/responsevector is independent of that feature. 
#'
#' @param datamatrix the data matrix with a column for each feature.
#' @param responsevector the vector of outcomes/responses (e.g. class labels). the length of this vector should match the amount of rows in datamatrix.
#' @param p.adj the adjustment method for the p-values. Any of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr' or 'none' are accepted.
#'
#' @return data with the features, their (adjusted) p-values and the linear estimate for the response, one for every column in the datamatrix.
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#' 
#' @examples
#' nSamples <- 10
#' nFeatures <- 20
#' data.matrix <- matrix( stats::runif(n=nFeatures*nSamples, min=0,max=100), 
#' ncol = nFeatures, nrow = nSamples)
#' 
#' response <- c( rep(0,nSamples/2), rep(1,nSamples/2) )
#' p_values <- relevant.features.p(datamatrix = data.matrix, responsevector = 
#' response, p.adj = 'none')
#' p_values_adjusted <- relevant.features.p( datamatrix = data.matrix, 
#' responsevector = response, p.adj = 'bonferroni')
#'
#' @export
#' 
#' @importFrom stats lm p.adjust var
#' 
relevant.features.p.extended <- function(datamatrix, responsevector, p.adj = "bonferroni") {
    
    ##### error checks #####
    if (!p.adj %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) 
        stop("p-value correction method is not one of the allowed possibilities")
    
    if (!class(responsevector) %in% "numeric") {
        print("responsevector was transformed to numeric")
        levels(responsevector) <- 1:length(levels(responsevector))
        responsevector <- as.numeric(responsevector)
    }
    
    if (length(responsevector) != nrow(datamatrix)) {
        stop("the length of the responsevector does not match with the amount of rows in the datamatrix")
    }
    ##### the actual function #####
    nfeat <- ncol(datamatrix)
    
    p <- rep(NA, nfeat)
    estimate.lm <- rep(NA, nfeat)
    for (coln in 1:nfeat) {
        if (stats::var(datamatrix[, coln]) == 0) {
            p[coln] <- 1
            estimate.lm[coln] <- 0 
        } else {
            linreg <- stats::lm(responsevector ~ datamatrix[, coln])
            test.lm <- summary(linreg)
            estimate.lm[coln] <- test.lm$coefficients[2,1]
            p[coln] <- summary(linreg)$coefficients["datamatrix[, coln]", "Pr(>|t|)"]
        }
    }
    
    p.adjusted <- stats::p.adjust(p, method = p.adj, n = nfeat)
    if (!is.null(colnames(datamatrix))) {
        features.indexes <- colnames(datamatrix)
    } else {
        features.indexes <- seq(1, length(p.adjusted))
    }
    
    P.data <- data.frame(matrix(NA, ncol = 2, nrow = length(p.adjusted)))
    colnames(P.data) <- c("index", "p.values")
    P.data$index <- features.indexes
    P.data$p.values <- p.adjusted
    P.data$estimates <- estimate.lm
    
    
    return(P.data)
    
}
