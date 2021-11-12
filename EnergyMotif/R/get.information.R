#' Information on Processes
#'
#' Function to extract the "normal" shape from the motifs discovered
#'
#' @param d original data frame, or data frame with at least a DateTime column
#' @param dmin list of not normalized data subsequences
#' @param starts starting points of the sequences in the original data
#' @param motifs detetcted motifs
#' @param max_length length to which the data has been normalized
#' @return returns a data frame with the information necessary for scheduling
#' @export

get.information <- function(d, dmin, starts, motifs, max_length) {
  # save the length of the processes
  length_processes <- unlist(lapply(dmin, function(x) length(x)))

  # get data frame for length information
  infos <- data.frame(Starts = starts)
  infos[["StartTime"]] <- d[["DateTime"]][starts]
  infos[["Length"]] <- length_processes
  infos[["Motif"]] <- NA
  infos[["Energy"]] <- unlist(lapply(dmin, sum))

  for (i in 1:length(motifs$Indices)) {
    infos[["Motif"]][(motifs$Indices[[i]] - 1)/max_length] <- i
  }
  return(infos)
}
