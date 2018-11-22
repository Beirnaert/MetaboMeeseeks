#' Helper function for geom_cobweb
#' 
#' This ggproto function draws the lines.
#'
#'
geom_cobweb_proto <- ggplot2::ggproto("geom_cobweb_proto", Geom,
required_aes = c("x", "y"),
default_aes = aes(colour = "black", alpha = 1, linetype = "solid", size = 2),
draw_key = draw_key_abline,
draw_group = function(data, panel_params, coord) {
    
    coords <- coord$transform(data, panel_params)
    
    ## Construct grid polygon
    grid::linesGrob(
        x= coords$x,
        y = coords$y,
        gp = grid::gpar(col = coords$colour,  
                        alpha = coords$alpha, 
                        lty = coords$linetype,
                        lwd = coords$size)
    )
    
}

)