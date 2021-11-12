#' SAX plot
#'
#' Visualize SAX of a time series.
#'
#' @param d a data frame with one column the time series to be plotted,
#' or a vector with the time series
#' @param x name of the column which includes the data, or new name
#' @param time specifiy the name of the column with the time data or set to "None"
#' @param a alphabet length for the SAX
#' @param path path where the plot should be saved
#' @param name name of the plot, default is SAX
#' @param normal specifiy the name of the column with the norm data
#' @param word complete word
#' @param pdf should there be a pdf file of the plot
#' @return returns the SAX plot
#' @export


vis.SAX <- function(d, x, a,
                    time = "None",
                    path,
                    name = "SAX",
                    normal,
                    word,
                    pdf = TRUE) {

  # create the break points form a normal distribution
  breaks <- round(stats::qnorm(p = seq(from = 0, to = 1, length.out = a + 1)),
                  digits = 2)

  if (time == "None") { # add a time column if none is existent
    d <- as.data.frame(d[ , c(x, normal)])
    names(d) <- c("Load", "Normal")
    d[["Time"]] <- 1:nrow(d)
  } else {
    d <- d[ , c(x, time, normal)]
    names(d) <- c("Load", "Time", "Normal")
  }

  # plot a step function
  p <- ggplot2::ggplot(d, ggplot2::aes(x = Time, y = Load)) +
    ggplot2::scale_x_continuous(expand = c(0, 0),
                                limits = c(-0.1, nrow(d) + 15)) +
    ggplot2::geom_line(data = d, ggplot2::aes(x = Time, y = Normal),
                       color = "gray70") +
    ggplot2::geom_step() +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black"),
                   axis.text = ggplot2::element_text(size = 18),
                   axis.title = ggplot2::element_text(size = 18, face = "bold"),
                   plot.title = ggplot2::element_text(size = 18, face = "bold"))

  # add lines for the alphabet to the step function plot
  alph <- letters[1:(length(breaks) - 1)]
  p <- p + ggplot2::geom_text(ggplot2::aes(nrow(d) + 5, -1.5 ,label = c("a"),
                                           vjust = -1), size = 8)

  # add lines and labels to lines
  for (i in 2:length(alph)) {
    p <- p + ggplot2::geom_hline(yintercept = breaks[i], linetype = "dotted") +
      ggplot2::geom_text(x = nrow(d) + 5,
                y = breaks[i],
                label = alph[i],
                vjust = -1,
                size = 8)
  }

  p <- p + ggplot2::ggtitle(paste0(word, collapse = ""))

  if (pdf == TRUE) {
    grDevices::pdf(file = paste0(path, name, ".pdf"),
                   width = 14/2.54, height = 10/2.54)

    print(p)

    graphics.off()
  } else {
    print(p)
  }
}
