# Script to extract the "interesting" sequences 
# MIT License
# 
# Copyright (c) 2021 Nicole Ludwig

# load the self-written EnergyMotif R-Package which includes most of the core
# functionality needed for eSAX
library(devtools)
library(roxygen2)
install("EnergyMotif")
library(EnergyMotif)

library(foreach)
library(ggplot2)
library(gtools)

data = NULL
peak = FALSE # option to calculate the peak power of the sequences

## load data and plot it -------------------------------------------------------
aggregation <- c("1 second", "2 min", "10 min", "2 hours")

for (a in aggregation) {
  load(paste0("data/", a, "/data_", a, ".RData"))

  names(data)[1:2] <- c("Time", "Power")
  x <- "Power"

  "%nin%" <- Negate("%in%")

  data[["Time"]] <- as.POSIXct(data[["Time"]])

  # plot the complete time series ----------------------------------------------

  p <- ggplot(data, aes(x = Time, y = Power)) +
    theme_bw()  +
    geom_step() +
    theme(legend.position = "none",
          panel.border    = element_blank(),
          panel.grid      = element_blank(),
          axis.title      = element_blank(),
          axis.line       = element_line(colour = "black"),
          axis.text       = element_text(size = 18))

  pdf.plot(p, n = paste0("CompleteTimeSeries_", a), w = 10, h = 8)

  # create the subsequences with the day or subday patterns --------------------
  # calculate the ecdf for the alphabet
  my.ecdf <- calc.ecdf(x    = data$Power,
                       plot = TRUE,
                       n    = paste0("ecdf_Power_", a, ".pdf"))

  save(my.ecdf, file = paste0("data/ecdf_", a, ".RData"))

  # create daily sequences
  window <- ifelse(a == "1 second", 24*60*60,
                   ifelse(a == "2 min", 24*60,
                          ifelse(a == "10 min", 24*(60/10),
                                 12)))
  # get the sequences
  dmin <- minimum.search(d      = data,
                         x      = "Power",
                         event  = "none",
                         window = window)

  # store the startpoints and sequences separately to not have lists of lists
  starts <- dmin[["startpoints"]]
  dmin   <- dmin[["dmin"]]

  # plot the sequences
  dat <- lapply(dmin, function(x) cbind(x = seq_along(x), y = x))
  # set names for the individual sequences
  list.names <- paste("Seq", as.character(seq(1, length(dmin))))

  # store the length off the individual sequences
  lns <- sapply(dat, nrow)

  # create a data frame with one y variable (all sequences after each other)
  # and a x variable (time steps)
  # add the name of the sequence as third column
  # by repeating the name according to the length of the sequence
  dat          <- as.data.frame(do.call("rbind", dat))
  names(dat)   <- c("Timesteps", "Load")
  dat$Sequence <- rep(list.names, lns)

  dat <- within(dat,
                Sequence <- factor(Sequence,
                                   levels = mixedsort(unique(dat$Sequence))))

  p <- ggplot(dat, aes(x = Timesteps, y = Load, colour = Sequence)) +
    theme_bw()  +
    geom_step() +
    facet_wrap(~ Sequence) +
    theme(legend.position = "none")

  pdf.plot(p, n = paste0("All_Sequences_",a), w = 14, h = 10)

  ts.subs <- dmin

  # this is the most important file for the later eSAX transformation
  # and motif discovery
  save(ts.subs, file = paste0("data/tssubs_", a, ".RData"))

}


## get the peak power ----------------------------------------------------------

if (peak == "TRUE") {
  load("data/tssubs_day.RData")

  splittime <- list()

  for (i in 1:length(ts.subs)) {
    y <- ts.subs[[i]]

    plot(y, type = "l")
    abline(h = quantile(y, 0.75), col = "red")

    cross <- which(y > quantile(y, 0.75))
    cross <- cross[which(cross > 60000)][1]

    abline(v = cross, col = "blue")

    splittime[[i]] <- data$Time[cross]
  }

  meantime <- mean(as.POSIXct(unlist(splittime),
                              origin = '1970-01-01',
                              tz     = "CET"))

  splitting <- which(data$Time == as.character(meantime))

  peak.sequences <- list()
  for (i in 1:length(ts.subs)) {
    peak.sequences[[i]] <- ts.subs[[i]][splitting:length(ts.subs[[i]])]
  }

  ts.subs <- peak.sequences

  save(ts.subs, file = "data/tssubs_peak.RData")
}

print("Done")

rm(list = ls())