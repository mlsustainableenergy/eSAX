#' Minimum Search
#'
#' Function to find Minima in the time series that could indicate Process starts
#'
#' @param d data frame
#' @param x name of the column containing the data to be evaluated
#' @param event (none, zero, minimum, custom) subsequences are either determined by a
#' minimum search or thorugh the points where they are zero or another
#' specified value. If none is selected the subsequences are predefined by
#' the window length
#' @param custom.event the customized value for the event start
#' @param window.size indicates the window size for the minimum search
#' @param window length of the subsequences if no event search is taking place,
#' if motif discovery follows it should be same window length
#' @return a list with the starting points and the subsequences
#' @export

minimum.search <- function(
  d,
  x,
  event = c("none", "zero", "minimum", "custom"),
  custom.event = 0.06,
  window.size = 100,
  window = NA) {

  if (event == "minimum") {
    cat("Searching for minima ... \n")
    # initialize vector for minima
    localmin <- c()

    # loop that finds all minima occuring in each run
    w <- window.size

    # find the minima in the window (use the first one found if more than one)
    for (i in 1:(nrow(d)/w)) {
      k <- i*w
      j <- (i*w) - w
      window <- d[j:k, x]
      localmin[i] <- c(which(window == min(window)) + ((i - 1)*w))[1]
    }

    cat("Preparing list ... \n")

    dmin <- list()

    dmin[[1]] <- d[1:localmin[1], x]

    for (i in 1:(length(localmin) - 1)) {
      dmin[[i + 1]] <- d[localmin[i]:(localmin[i + 1] - 1), x]
    }

    dmin[[length(dmin) + 1]] <- d[localmin[length(localmin)]:nrow(d), x]

  } else if (event == "zero") {
    cat("Searching for zeros ... \n")

    localmin <- c()
    start <- which(d[[x]] == 0)

    for (i in 1:(length(start) - 1)) {
      if (d[(start[i] + 1), x] != 0) {
        localmin[[length(localmin) + 1]] <- start[i]
      }
    }

    endpoints <- c()

    # find next point where it is zero again
    for (i in 1:length(localmin)) {
      endpoints[[length(endpoints) + 1]] <- start[which(start > localmin[i])[1]]
    }

    localmin <- c(localmin, endpoints)
    localmin <- localmin[order(stats::na.omit(localmin))]

    dist.to.min <- c()

    for (i in 2:length(localmin)) {
      dist.to.min[i - 1] <- localmin[i] - localmin[i - 1]
    }

    cat("Preparing list ... \n")

    dmin <- list()
    for (i in 1:(length(localmin) - 1)) {
      dmin[[i]] <- d[localmin[i]:localmin[i + 1], x]
    }
    a <- c()
    for (i in 1:length(dmin)) {
      if (length(which(dmin[[i]] != 0)) == 0) {
        a[[length(a) + 1]] <- i
      }
    }
    dmin <- dmin[-a]

  } else if (event == "custom") {
    cat("Searching for custom event ... \n")

    localmin <- c()
    start <- which(d[[x]] <= custom.event)

    for (i in 1:(length(start) - 1)) {
      if (d[(start[i] + 1), x] != custom.event) {
        localmin[[length(localmin) + 1]] <- start[i]
      }
    }

    endpoints <- c()

    # find next point where it is zero again
    for (i in 1:length(localmin)) {
      endpoints[[length(endpoints) + 1]] <- start[which(start > localmin[i])[1]]
    }

    localmin <- c(localmin, endpoints)
    localmin <- localmin[order(stats::na.omit(localmin))]

    dist.to.min <- c()

    for (i in 2:length(localmin)) {
      dist.to.min[i - 1] <- localmin[i] - localmin[i - 1]
    }

    cat("Preparing list ... \n")

    dmin <- list()
    for (i in 1:(length(localmin) - 1)) {
      dmin[[i]] <- d[localmin[i]:localmin[i + 1], x]
    }
    a <- c()
    for (i in 1:length(dmin)) {
      if (length(which(dmin[[i]] != 0)) == 0) {
        a[[length(a) + 1]] <- i
      }
    }
    if (length(a) > 0) {
      dmin <- dmin[-a]
    }


  } else if (event == "none") {
    cat("Preparing subsequences ... \n")

    dmin <- list()
    # store the subsequences of size window length for Motif discovery in dmin
    for (i in 1:floor(nrow(d)/window)) {
      dmin[[i]] <- d[((i - 1)*window + 1):(i*window), x]
    }

    localmin <- c()
    # save the startpoints (window length distance)
    for (i in 1:length(dmin)) {
      localmin[i] <- (i - 1)*window + 1
    }

    cat("Preparing list ... \n")

  }

  return(list(dmin = dmin, startpoints = localmin))
}
