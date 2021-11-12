#' Uniform
#'
#' Function to transform data to uniform data
#'
#' @param x original vector
#' @return vector with uniform data
#' @export

uniform <- function(x){
  (x - min(x))/(max(x) - min(x))
}