#' Get the process information
#'
#' Wrapper function for complete motif discovery for energy time series data
#'
#' @param d original data
#' @param x name of the column which includes the data
#' @param plot.var define which variable should be plotted
#' @param machine name of the machine which is considered
#' @param name name of the time series under investigation
#' @param event which event search method is used for the subsequences ("none",
#' zero", "minimum" or "custom")
#' @param custom.event custom event start indicator
#' @param norm normalise the subsequences to a common length
#' @param window.size window size for minimum search
#' @param alphabet define the alphabet size for the motif discovery
#' @param w word size for motif discovery
#' @param window window size for motif discovery
#' @param path path to where the plots are printed and the data is saved
#' @param cr count ration for motif discovery (default is 1.5)
#' @param mr value used to add other possible members to a motif candidate
#' @param mask mask size used for random projection
#' @param epsilon epsilom used for threshold in motif discovery
#' @param iter iterations used in random projection algorithm
#' @return returns complete process motif information and vector with result
#' @export

GetProcesses <- function(d,
                         x = "median",
                         plot.var = "Median",
                         machine,
                         name,
                         event = c("none", "zero", "minimum", "custom"),
                         custom.event = 0,
                         norm = TRUE,
                         window.size = 100,
                         alphabet,
                         w = NA,
                         window = NA,
                         cr = 1.5,
                         mr = 1.2,
                         mask = 2,
                         iter = 25,
                         epsilon = 0.1,
                         path) {

  minima <- minimum.search(d, x, event, window.size, window = window)

  dmin <- minima[["dmin"]]

  if (is.data.frame(dmin[[1]])) {
    dmin <- lapply(dmin, function(x) as.numeric(unlist(x)))
  }

  # save the start points
  starts <- minima[["startpoints"]]

  if (event != "none") { # clean dmin if event search was done

    # only keep those sequences which do not consist of only zero entries
    keep <- which(lapply(dmin,
      function(x) (length(which(x == 0)) == length(x))) == FALSE)

    dmin <- dmin[keep]
    starts <- starts[keep]

    # delete the subsequences which are shorter than 5
    too.short <- which(lapply(dmin, function(x) length(x)) <= 5)

    # need if here in case there are no short sequences
    if (length(too.short) > 0) {
      dmin <- dmin[-too.short]
      starts <- starts[-too.short]
    }

    # normalize the data subsequences to a common length
    if (norm == TRUE) {
      cat("Normalize Data \n")
      samelength <- normalize.length(dmin, length = "Quantile",
        q = 0.9,
        plot = TRUE,
        graphicspath = path)

      max_length <- length(samelength[[1]])
      newtimeseries <- unlist(samelength)

      # assign the word length according to max_length if not defined before
      w <- ifelse(is.na(w), max_length, w)

      # same for window length
      window <- ifelse(is.na(window), max_length, window)

      # save the new time series as R data frame and csv
      save(newtimeseries, file = paste0(path, machine, "-",
        max_length, "_ts.RData"))

      utils::write.csv(newtimeseries,
        file = paste0(path, machine, "-", max_length, "_ts.csv"))

    } else {# if no normalisation should be done on the time series, keep original

      newtimeseries <- unlist(dmin)

    }
  } else if (event == "none") {

    newtimeseries <- d[, x]
  }

  # initialise the top paramter (needed for stopping criteria later)
  top <- 1

  # assign count ratio (stopping criteria)
  #assign("my.cr", cr, envir = globalenv())

  #repeat {
    tryCatch({
      cat("searching Motif with count ratio", cr, "and alphabet size",
          alphabet, "\n")

      allmotifs <-  motif.discover(newtimeseries,
                                   global.norm = TRUE,
                                   local.norm = TRUE,
                                   window.size = window,
                                   overlap = 0,
                                   w = w,
                                   a = alphabet,
                                   mask.size = mask,
                                   eps = epsilon,
                                   iter = iter,
                                   count.ratio.1 = cr,
                                   max.dist.ratio = mr)

      top <- length(unlist(allmotifs$Indices))

    #if (top != 1) {
      #break
    #}
  #}

  # plot the times series
  motif.plot(x = newtimeseries, w = w, window = window,
             motif = allmotifs,
             graphicspath = paste0(path, top, "_P", "_"),
             pdf = TRUE)

  save(allmotifs, file = paste0(path, top, "_new_found_motifs.RData"))

  if (event == "none") {
    max_length <- window
  }

  infos <- get.information(d, dmin, starts, motifs = allmotifs, max_length)

  utils::write.csv(infos, file = paste0(path, top, "_Information_", name,
                                        ".csv"))

  cat("Getting standard shape for Motifs ... \n")

  get.shape(motif = allmotifs,
            ts = newtimeseries,
            plot = TRUE,
            save = TRUE,
            all = TRUE,
            infos = infos,
            plot.var = plot.var,
            machine = machine,
            top = top,
            path = path)

  for (i in 1:length(allmotifs$Indices)) {
    load(paste0(path, top, "_aa_", i, "df.Rdata"))
    aa <- aa[,c("Q20", "Mean", "Median", "Q80")]
    utils::write.csv(aa, file = paste0(top, "_Process_", i, ".csv"))
  }

  result <- c("name" = name,
              "window" = window,
              "alphabet" = alphabet,
              "words" = w,
              "count.ratio 1" = cr,
              "max to dist ratio" = mr,
              "sequence length" = length(allmotifs$Indices),
              "# motifs" = top,
              "mask size" = mask,
              "iter" = iter,
              "eps" = epsilon)

    }, error = function(e){
      cat("No Motifs found with selected parameters", "\n")

      result <- c("name" = name,
                  "window" = window,
                  "alphabet" = alphabet,
                  "words" = w,
                  "count.ratio 1" = cr,
                  "max to dist ratio" = mr,
                  "sequence length" = NA,
                  "# motifs" = NA,
                  "mask size" = mask,
                  "iter" = iter,
                  "eps" = epsilon)

      write.table(result, file = paste0("Error_", i, ".txt"))
      #assign("my.cr", my.cr + 0.5, envir = globalenv())
    })

  cat("Done \n")

  return(result)
}