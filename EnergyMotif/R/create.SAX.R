#' Perfom a SAX for time series
#'
#' The function create SAX symbols for a univariate time series.
#' The details of this method can be referred to J. Lin, E. Keogh, L. Wei, S. Lonardi.
#' Experiencing SAX: a novel symbolic representation of time series
#' Code as in the package TSMining but with a longer alphabet
#'
#' @param x a numeric vector representing the univarate time series
#' @param w defines the word size used for SAX transformation
#' @param a defines the alphabet size used for SAX transformation
#' @param eps s the minimum threshold for variance in subsequence and should be
#' a numeric value. If the subsequence considered has a smaller variance than eps,
#' it will be represented as a word using the middle alphabet. The default value is 0.1
#' @param norm logical, deciding whether standardization should be applied to x.
#' If True, x is standardized using mean and standard deviation
#' @param return.pieces if TRUE also returns the pieces of the PAA
#' @return The function returns a SAX representation of x.
#' @export


create.SAX  <- function(x, w, a, eps, norm, return.pieces = FALSE) {

  require(foreach)

  alphabet <- c(letters, paste0(letters,letters), paste0(letters,letters, letters))

  i = NULL
  if (stats::sd(x) <= eps) {
    sym <- rep(alphabet[round((1 + a)/2, digits = 0)], w)
    } else {
      # Normalize the data to have 0 mean and 1 standard deviation
      # before piecewise aggregation
    if (norm == TRUE) {
      data.nor <- (x - mean(x))/stats::sd(x)
      } else {
      data.nor <- x
    }

    # Perform the piecewise aggregation
    ind <- round(seq(from = 1,
                     to = length(data.nor), length.out = w + 1),
                 digits = 0)
    pieces <- foreach::foreach(i = 1:(length(ind) - 1), .combine = c) %do% {
      if (i != (length(ind) - 1)) {
        piece <- data.nor[ind[i]:(ind[i + 1] - 1)]
      } else {
          piece <- data.nor[ind[i]:ind[i + 1]]
          }
      return(mean(piece, na.rm = T))
    }

    #Perform alphabet assignment
    let <- alphabet[1:a]
    #Create breaks points based on Gaussian normal distribution
    bks <- round(stats::qnorm(p = seq(from = 0, to = 1, length.out = a + 1)),
                 digits = 2)
    sym <- foreach::foreach(i = 1:length(pieces), .combine = c) %do% {
      obs <- pieces[i]
      let[max(which(bks < obs))]
    }
  }

  if (return.pieces == FALSE) {
    return(sym)
  } else {
    return(list(sym = sym, pieces = pieces, ind = ind, bks = bks))
  }

}
