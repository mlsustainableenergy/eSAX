#' Function which uses a stored plot and puts it in a pdf
#'
#' @name pdf.plot
#' @param p the graph to be saved in a pdf
#' @param n name for the pdf file
#' @param w width of the pdf file
#' @param h height of the pdf file
#' @return pdf file of plot
#' @export

pdf.plot <- function(p, n = "graph", w = 14/2.54, h = 10/2.54){

  grDevices::pdf(file = paste0(n, ".pdf"),
                 width = w,
                 height = h)

  plot(p)

  graphics.off()
}