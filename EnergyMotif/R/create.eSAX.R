#' Perfom an eSAX transformation for a time series
#'
#' The function creates eSAX symbols for an univariate time series.
#'
#' @param x a numeric vector representing the univarate time series
#' @param w defines the word size used for eSAX transformation
#' @param b the breakpoints used for the eSAX representation
#' @return The function returns an eSAX representation of x.
#' @export

create.eSAX <- function(x, b, w, ecdf) {

  i = NULL
  # Perform the piecewise aggregation
  ind <- round(seq(from = 1, to = length(x), length.out = w + 1), digits = 0)

  pieces <- foreach::foreach(i = 1:(length(ind) - 1), .combine = c) %do% {
    if (i != (length(ind) - 1)) {
      piece <- x[ind[i]:(ind[i + 1] - 1)]
    } else {
      piece <- x[ind[i]:ind[i + 1]]
    }
    return(mean(piece, na.rm = T))
  }

  # create an alphabet with douple and triple letter combinations (a, aa, aaa)
  alphabet <- c(letters, paste0(letters,letters), paste0(letters,letters, letters))

  # assign the alphabet
  let <- alphabet[1:length(b)]

  # add symbols to sequence according to breakpoints
  sym <- foreach::foreach(i = 1:length(pieces), .combine = c) %do% {
    obs <- pieces[i]
    let[max(which(b <= obs))]
  }

  return(list(sym = sym, pieces = pieces, ind = ind, bks = b))
}