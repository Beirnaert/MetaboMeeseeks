#' Cobweb plot
#' 
#' The cobweb geom is an addition to the \pkg{ggplot2} package for plotting
#' cobweb or spider plots.
#' 
#' @param mapping (from ggplot call) A mapping of ggplot variables. geom_cobweb
#'   expects an AUCs variable.
#' @param data (from ggplot call) The data frame supplied to the
#' @param stat (from ggplot call) ggplot identity stat.
#' @param position (from ggplot call) ggplot identity position.
#' @param na.rm (from ggplot call) remove NAs, default is FALSE.
#' @param theme Default theme settings for correct visualization.
#' @param show.legend Whether to show the legend. Default is TRUE.
#' @param grid_lty The linetype of the background grid. Deafult is 'dashed'.
#' @param grid_lwd The linewidth of the background grid. Deafult is 1.
#' @param grid_alpha The alpha of the background grid. Deafult is 0.75.
#' @param grid_fontsize The fontsize of the grid labels. Deafult is 12.
#' @param grid_Nticks The number of gridlines (makes use of
#'   \code{\link[base]{pretty}} so number is not a hard threshold). Default is
#'   4.
#' @param grid_ticks vector with values of gridlines to be drawn (unlike
#'   grid_Nticks, this is exact).
#' @param grid_labels The labels to use for the gridlines.
#' @param grid_label_position The position alonng the cobweb where to plot the
#'   grid labels. (numeric value)
#' @param grid_label_hjust The hjust of grid_labels
#' @param grid_label_vjust The vjust of grid_labels
#' @param ... Extra aesthetics parameters, see Aesthetics section.
#'   
#' @section Aesthetics: \code{geom_cobweb} understands the following aesthetics
#'   (required aesthetics are in bold): 
#'   \itemize{ 
#'   \item \strong{\code{AUCs}} 
#'   \item \code{alpha} 
#'   \item \code{colour} 
#'   \item \code{group} 
#'   \item \code{linetype} 
#'   \item \code{size} }
#'   
#'   
#' @export
#' @examples
#' library(ggplot2)
#' n_comparison <- 3
#' AUC.df <- data.frame(type = c(rep("ROC",6),rep("Random",6)),  
#'                     AUCs = c(0.34, 1.00, 0.56, 0.40, 0.37, 0.45, rep(1/n_comparison, 6)))
#' 
#' ggplot(AUC.df, aes(AUCs = AUCs,  colour = type)) + 
#'     geom_cobweb() + 
#'     theme_cobweb() +
#'     ggtitle("Testen") +
#'     coord_equal()
#'  
#'  
#'  
#' ggplot(AUC.df, aes(AUCs = AUCs,  colour = type)) + 
#'     geom_cobweb(grid_Nticks = 2,
#'                 show.legend = TRUE,
#'                 grid_label_position = 2) + 
#'     theme_cobweb() +
#'     ggtitle("Testen") +
#'     coord_equal()
#'  
geom_cobweb <- function(mapping = NULL, data = NULL, stat = 'identity',
                         position = 'identity', na.rm = FALSE,
                         show.legend = TRUE, 
                         theme = list(panel.background = element_rect(fill = "white", colour = "white")), 
                         grid_lty = "dashed", 
                         grid_lwd = 1, 
                         grid_alpha = 0.75, 
                         grid_fontsize = 12,
                         grid_Nticks = 4,
                         grid_ticks = NULL, 
                         grid_labels = NULL, 
                         grid_label_position = 1,
                         grid_label_hjust = 0.5,
                         grid_label_vjust = -0.5, ...) {
    
    
    convert_to_cobweb_bg <- function(data){
        if(!"AUCs" %in% colnames(data)){
            stop("Haha such silly code, it requires 'AUCs' to be present in the variable names (for the values to plot)")
        }
        if(ncol(data) == 2){
            N <- unique(summary(data[,which(colnames(data) != "AUCs")[1]]))
        }else if(ncol(data) == 1){
            N <- nrow(data)
        }else{
            factor_vars <- NULL
            col_vars <- lapply(data,class)
            for(k in 1:length(col_vars)){
                if("factor" %in% col_vars[k]){
                    factor_vars <- c(factor_vars, k)
                }
            }
            N <- unique(summary(data[,factor_vars[1]]))
            message(paste(colnames(data)[factor_vars[1]], "chosen as class variable", sep = " "))
        }
        
        maxAUC <- max(data$AUCs)
        if(is.null(grid_ticks)){
            ticks <- pretty(c(0,maxAUC), n = grid_Nticks)
            labels <- as.character(ticks)
        } else if(!is.null(grid_ticks) & is.null(grid_labels)){
            ticks <- as.numeric(grid_ticks)
            labels <- as.character(grid_ticks)
        } else if(!is.null(grid_ticks) & !is.null(grid_labels)){
            ticks <- as.numeric(grid_ticks)
            if(length(grid_labels) != length(grid_ticks)){
                labels <- as.character(grid_ticks)
            } else{
                labels <- as.character(grid_labels)
            }
        }
        
        if(0 %in% ticks){
            labels <- labels[ticks !=0]
            ticks <- ticks[ticks !=0]
        }
        
        background_data <- data.frame(ticks = ticks, labels = labels, N = N)
        return(background_data)
    }
    
    stat_cobweb_bg <- ggproto("stat_cobweb_bg", Stat,
                              compute_group = function(data, scales) {
                                  
                                  ticks <- rep(data$ticks, data$N)
                                  
                                  bg_data <- data.frame(x = rad2polar2xy(ticks)$x, 
                                                        y = rad2polar2xy(ticks)$y, 
                                                        labels = data$labels[1])
                                  bg_data
                              },
                              required_aes = c("ticks", "labels", "N")
    )
    
    stat_cobweb <- ggproto("stat_cobweb", Stat,
                           compute_group = function(data, scales) {
                               
                               data_df <- data.frame(x = rad2polar2xy(data$AUCs)$x,
                                                     y = rad2polar2xy(data$AUCs)$y)
                               data_df
                           },
                           required_aes = c("AUCs")
    )
    
    
    
    list(
        ggplot2::layer(
            geom = geom_cobweb_bg_proto, # step 3 (gets a layer)
            mapping = aes_string(ticks = "ticks", labels = "labels", N = "N"),#aes(ticks = ticks, labels = labels), # goes into geom (step 3)
            data = convert_to_cobweb_bg, # step 1
            stat = stat_cobweb_bg, # step 2 (gets a layer) so this should get 1 line and replicate that N times
            position = position,
            show.legend = FALSE,
            inherit.aes = FALSE,
            params = list(na.rm = na.rm,
                          linetype = grid_lty,
                          size = grid_lwd,
                          alpha = grid_alpha,
                          fontsize = grid_fontsize,
                          label_position = grid_label_position,
                          label_hjust  = grid_label_hjust,
                          label_vjust = grid_label_vjust)
        ),
        
        ggplot2::layer(
            geom = geom_cobweb_proto, 
            mapping = mapping,
            data = data, 
            stat = stat_cobweb, 
            position = position, 
            show.legend = show.legend, 
            inherit.aes = TRUE,
            params = list(na.rm = na.rm,...)
        ) 
    )
    
}
