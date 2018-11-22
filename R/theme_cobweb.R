#' Theme for cobweb visualization
#'
#' \code{theme_cobweb} Provides the theme to the geom_cobweb ggplot2 function. It is advised
#' to also add coord_equal() to the plot to achieve a cobweb with of equal lengths. 
#' 
#'  
#' @param base_size base font size
#' @param base_family base font family
#'
#' @examples
#' library(ggplot2)
#' n_comparison = 3
#' AUC.df = data.frame(type = c(rep("ROC",6),rep("Random",6)),  
#'                     AUCs = c(0.34, 1.00, 0.56, 0.40, 0.37, 0.45, rep(1/n_comparison, 6)))
#' ggplot(AUC.df, aes(AUCs = AUCs,  colour = type)) +
#'     geom_cobweb() + 
#'     theme_cobweb() +
#'     ggtitle("Test plot") 
#'
#' @name theme_cobweb
#' 
#' @import ggplot2
#' 
#' @export
theme_cobweb <- function (base_size = 11, base_family = "")
{
    theme_grey(base_size = base_size, base_family = base_family) %+replace%
        theme(
            plot.background = element_rect(fill = "white", colour = "white"),
            panel.background = element_rect(fill = "white", colour = "white"),
            #legend.box.spacing = grid::unit(-0.6, unit="cm"),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            
            legend.key = element_rect(fill = NA, colour = "white"),
            
            panel.border = element_blank(),
            panel.grid = element_blank(),
            
            plot.title = element_text(hjust = 0.5,
                                      size = base_size+7)
        ) 
}
