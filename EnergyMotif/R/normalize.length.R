#' Normalize length of time series
#'
#' Function to find Minima in the time series that could indicate Process starts
#'
#' @param dmin list with time series of different length
#' @param length length of the normalized time series, "Max" for maximum length
#' "Quantile" for a quantile which should be used, "Value" specify a length
#' @param q if length is quantile, then we need q in [0,1]
#' @param value of length method is value, need to specify here
#' @param plot TRUE/FALSE, should there be a plot of all sequences
#' @param graphicspath path to where the normalized time series plot is saved
#' @param list (TRUE/FALSE), not sure what I wanted to do there
#' @return list of the length normalized time series
#' @export

normalize.length <- function(dmin, length = c("Max", "Quantile", "Value"),
                             q = 0.9,
                             value = 100,
                             plot = TRUE,
                             graphicspath = graphicspath,
                             list = TRUE) {

  if (plot == TRUE) {
    grDevices::pdf(file = paste0(graphicspath, "Not-Normalized.pdf"),
                   height = 8,
                   width = 12,
                   onefile = TRUE)

    for (i in 1:length(dmin)) {

      plot(dmin[[i]],
           type = "l",
           xlab = "Time", ylab = "Load",
           main = "")
    }

    grDevices::graphics.off()
  }

  # save the length of the processes
  length_processes <- unlist(lapply(dmin, function(x) length(x)))

  # normalized to which length
  if (length == "Max") {
    max_length <- max(length_processes)
  } else if (length == "Quantile") {
    max_length <- stats::quantile(length_processes, q)
    max_length <- sort(length_processes[which(length_processes > max_length)])[1]
  } else {
    max_length <- value
  }

  samelength <- list()
  p.maxlength <- which(length_processes == max_length)
  samelength[[1]] <- unlist(dmin[p.maxlength[1]])

  p <- which(length_processes != max_length)

  output <- 1:max_length

  for (j in p) {
    input <- unlist(dmin[[j]])
    input <- 1:length(input)
    interp_output <- stats::approx(x = seq(0, max(output), length.out = length(input)),
                            y = input,
                            xout = output)
    newtimes <- round(interp_output$y)
    newtimes[1] <- 1

    samelength[[length(samelength) + 1]] <- dmin[[j]][newtimes]
  }

  if (plot == TRUE) {
    grDevices::pdf(file = paste0(graphicspath, "Normalized.pdf"),
        height = 8,
        width = 12,
        onefile = TRUE)

    for (i in 1:length(samelength)) {

      plot(samelength[[i]],
           type = "l",
           xlab = "Time", ylab = "Load",
           main = "")
    }

    grDevices::graphics.off()
  }

  return(samelength)
}

