% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.shape.R
\name{get.shape}
\alias{get.shape}
\title{Get the standard shape of each motif found in time series}
\usage{
get.shape(motif, ts, plot = TRUE, save = TRUE, all = TRUE, noise = TRUE,
  infos = infos, plot.var = "Median", machine = machine, top, path,
  publication = TRUE)
}
\arguments{
\item{motif}{output of motif discovery}

\item{ts}{time series which was used for motif discovery}

\item{plot}{(logical) should there be a plot}

\item{save}{(logical) should the data be saved}

\item{all}{(logical) should a single plot with all motifs be plotted}

\item{noise}{(logical) should a noise motif be plotted (needs info file at the moment)}

\item{infos}{data frame with information about which motif starts at which
time step}

\item{plot.var}{which variable is plotted comparing all motifs}

\item{machine}{name of current machine under investigation}

\item{top}{count of motifs found}

\item{path}{path to where the plot is printed}

\item{publication}{changes font size for graphics}
}
\value{
returns a list with data frames of "normal" shapes with mean,
median, q20 and q80
}
\description{
Function to extract the "normal" shape from the motifs discovered
}
