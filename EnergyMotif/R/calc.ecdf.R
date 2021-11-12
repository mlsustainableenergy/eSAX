#' Calculated the empirical CDF of a time series
#'
#' @param x a numeric vector representing the univarate time series
#' @param plot indicates if the ecdf should be plotted
#' @param n name of the plot
#' @return The function returns a ecdf function for the time series
#' @export

calc.ecdf <- function(x, plot = TRUE, n = 'ecdf.pdf') {
  # only use the values which are not zero for the calculation of the ecdf
  non.null <- x[which(x != 0)]

  # calculate the ecdf
  my.ecdf <- stats::ecdf(non.null)

  if (plot == TRUE) {
    # plotting the empirical cumulative distribution function
    # print a pdf image to a desired folder
    grDevices::pdf(n)

    graphics::plot(my.ecdf, xlab = 'Sample Quantiles',
                   ylab = '',
                   main = 'Empirical Cumluative Distribution')

    # add label for y-axis
    # the "line" option is used to set the position of the label
    # the "side" option specifies the left side
    graphics::mtext(text = expression(hat(F)[n](x)), side = 2, line = 2.5)
    grDevices::dev.off()
  }
  return(my.ecdf)
}