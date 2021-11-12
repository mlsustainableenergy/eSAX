#' Nice plot of the motifs discovered
#'
#' This functions creates a nice ggplot of the disovered motifs highlighted in
#' the original time series and plotted individually (if wanted)
#'
#' @name motif.plot
#' @param x the original time series data frame
#' @param w word length used for motif discovery
#' @param window the window size used for the motif discovery
#' @param motif the result of the motif.discover function
#' @param pdf is plotted as pdf or in the viewer
#' @param graphicspath path where to put the pdf plot if necessary
#' @param datestart beginning of date to be show on x axis
#' @param dateformat format of the date for the x axis
#' @param dateby steps in the date function (default 1 hour)
#' @param subfigure.1 plot 1 is used is a subfigure, thus the font size is bigger
#' @param subfigure.2 plot 2 is used is a subfigure, thus the font size is bigger
#' @return nice ggplots
#' @export

motif.plot <- function(x,
                       w,
                       window,
                       motif,
                       pdf = TRUE,
                       graphicspath,
                       datestart = "01.01.2015 00:00:00",
                       dateformat = "%d.%m.%Y %H:%M:%S",
                       dateby = "1 hour",
                       subfigure.1 = FALSE,
                       subfigure.2 = TRUE) {

  vis <- TSMining::Func.visual.SingleMotif(single.ts = x,
                                           window.size = window,
                                           motif.indices = motif$Indices)

  vis[[1]] <- stats::na.omit(vis[[1]])

  for (i in 1:length(vis[[2]])) {
    vis[[2]][[i]] <- stats::na.omit(vis$data.2[[i]])
  }

  myDate <- timeDate::timeSequence(from = datestart,
                                   format = dateformat,
                                   by = dateby,
                                   length.out = nrow(vis$data.1))

  vis$data.1$date <- as.POSIXct(myDate)

  if (length(motif$Indices) > 1) {
    ind <- c(motif$Indices[[1]], motif$Indices[[2]])
  } else {
    ind <- motif$Indices[[1]]
  }

  ind2 <- c(ind + window)
  ind <- c(ind, ind2)

  x <- as.data.frame(x)
  names(x) <- "Load"
  x[["DateTime"]] <- as.POSIXct(myDate)

  if (subfigure.1 == FALSE) {
    s <- 24
  } else if (subfigure.1 == TRUE) {
    s <- 32
  }

  e <- ggplot2::ggplot(x, ggplot2::aes(DateTime, Load)) +
      ggplot2::geom_line() +
      ggplot2::labs(y = "Load in MW", x = "Time Steps") +
      ggplot2::theme_bw() +
      ggplot2::geom_vline(xintercept = as.numeric(x[ind,"DateTime"])) +
      ggplot2::theme(panel.border = ggplot2::element_blank(),
                     axis.title = ggplot2::element_text(size = s),
                     legend.text = ggplot2::element_text(size = s),
                     axis.text.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     legend.title = ggplot2::element_text(size = s))

  if (pdf == TRUE) {
    pdf(file = paste0(graphicspath, "Motifs_", window, "_", w, ".pdf"),
        height = 8,
        width = 12,
        colormodel = "gray")

    print(e)

    grDevices::dev.off()
  } else {
    print(e)
  }

  if (subfigure.2 == FALSE) {
    ss <- 24
  } else if (subfigure.2 == TRUE) {
    ss <- 32
  }

  if (length(vis[[2]]) > 0) {
      for (l in 1:length(vis[[2]])) {
        f <- ggplot2::ggplot(data = vis$data.2[[l]]) +
      ggplot2::geom_line(ggplot2::aes(x = Time, y = Value, linetype = Instance)) +
      ggplot2::theme_bw() +
      ggplot2::scale_colour_grey() +
      ggplot2::labs(y = "Load in MW", x = "Time steps") +
      ggplot2::theme(panel.border = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.title = ggplot2::element_text(size = ss),
            legend.position = "none"
      )

    if (pdf == TRUE) {
      pdf(file = paste0(graphicspath, "Motifs_", window, "_", w, "_", l, ".pdf"),
          height = 8,
          width = 12,
          colormodel = "gray")

      print(f)

      grDevices::dev.off()
    } else {
      print(f)
    }
      }
  }
}

