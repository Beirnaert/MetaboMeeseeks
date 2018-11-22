#' Helper function for geom_cobweb
#' 
#' This ggproto function draws the background stuff.
#'
#'
geom_cobweb_bg_proto <- ggplot2::ggproto("geom_cobweb_bg_proto", Geom,
                                         required_aes = c("x", "y"),
                                         default_aes = aes(colour = "black", alpha = 0.75, linetype = "solid", 
                                                           size = 1, fontsize = 12, label_position = 1,
                                                           label_hjust = 0.5, label_vjust = -0.5),
                                         draw_key = draw_key_abline,
                                         draw_group = function(data, panel_params, coord) {
                                             
                                             coords <- coord$transform(data, panel_params)
                                             
                                             ## Construct grid polygon
                                             bg_grob <- grid::linesGrob(
                                                 x = coords$x,
                                                 y = coords$y,
                                                 gp = grid::gpar(col = coords$colour,  
                                                                 alpha = coords$alpha, 
                                                                 lty = coords$linetype,
                                                                 lwd = coords$size)
                                             )
                                             
                                             lab_position <- length(coords$x) - (coords$label_position[1] %% (length(coords$x)-1))
                                             
                                             bg_labels <- grid::textGrob(label = as.character(data$labels[1]),
                                                                         x = coords$x[lab_position],
                                                                         y = coords$y[lab_position],
                                                                         hjust = data$label_hjust[1],
                                                                         vjust = data$label_vjust[1],
                                                                         default.units = "npc",
                                                                         gp=grid::gpar(alpha = coords$alpha[1], fontsize = coords$fontsize))
                                             
                                             grid::gList(bg_grob,bg_labels)
                                             
                                         }
                                         
)