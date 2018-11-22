#' Helper function for polar plotting
#' 
#' \code{rad2polar2xy} Converts single values to polar coordinates, with the points evenly spaced out.
#' These points are then transformed to cartesian coordinates for plotting with geom_cobweb.
#'
#' @param x The values to convert to evenly spaced polar coordinates, then to x,y coordinates.
#'
rad2polar2xy <- function(x){
    thetas <- seq(pi/2, 5*pi/2, length.out = length(x)+1 )
    radii <- c(x, x[1])
    
    coordinates <- data.frame(x = 0.5+(radii*cos(thetas)/2), y = 0.5+(radii*sin(thetas)/2))
    coordinates[nrow(coordinates),] <- coordinates[1,]
    return(coordinates)
}