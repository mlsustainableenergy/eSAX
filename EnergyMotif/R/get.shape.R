#' Get the standard shape of each motif found in time series
#'
#' Function to extract the "normal" shape from the motifs discovered
#'
#' @param motif output of motif discovery
#' @param ts time series which was used for motif discovery
#' @param plot (logical) should there be a plot
#' @param save (logical) should the data be saved
#' @param infos data frame with information about which motif starts at which
#' time step
#' @param noise (logical) should a noise motif be plotted (needs info file at the moment)
#' @param all (logical) should a single plot with all motifs be plotted
#' @param plot.var which variable is plotted comparing all motifs
#' @param machine name of current machine under investigation
#' @param top count of motifs found
#' @param path path to where the plot is printed
#' @param publication changes font size for graphics
#' @return returns a list with data frames of "normal" shapes with mean,
#' median, q20 and q80
#' @export


get.shape <- function(motif,
  ts,
  plot = TRUE,
  save = TRUE,
  all = TRUE,
  noise = TRUE,
  infos = infos,
  plot.var = "Median",
  machine = machine,
  top,
  path,
  publication = TRUE) {

  # transform list to data frame (current process length)
  mm <- motif
  motifdata <- as.data.frame(ts)
  names(motifdata) <- "Value"
  motifdata[["Time"]] <- seq(1:nrow(motifdata))

  l <- dim(mm$Subs)[2] - 2
  xx <- list()

  normaleshape <- list()

  for (k in 1:length(mm$Indices)) {
    s <- mm$Indices[[k]][1]

    # initialize extra data frame with first instance of motif
    xx[[k]] <- as.data.frame(motifdata[s:(s + l), "Value"])

    # add all other instances of the motif to the data frame
    for (y in 2:length(mm$Indices[[k]])) {
      b <- mm$Indices[[k]][y]
      xx[[k]] <- as.data.frame(cbind(xx[[k]],
        motifdata[b:(b + l), "Value"]))
    }
  }

  if (noise == TRUE) {
    # add noise motif
    Noise <- unlist(which(is.na(infos[["Motif"]])))
    Noise.List <- list()
    seq.length <- nrow(xx[[1]])
    for (i in 1:length(Noise)) {
      start <- (Noise[i] - 1)*seq.length + 1
      Noise.List[[i]] <- as.data.frame(ts[start:(start + (seq.length - 1))])
    }

    notgood <- length(which(
      unlist(lapply(Noise.List,
        function(x) length(which(is.na(x))))) == nrow(xx[[1]])))

    if (notgood != length(Noise.List)) {
      xx[[length(xx) + 1]] <- as.data.frame(Noise.List)
      no.noise = FALSE
    } else {
      no.noise <- TRUE
    }
  }


  # give all the instances a letter
  for (r in 1:length(xx)) {
    names(xx[[r]])[1:length(xx[[r]])] <- c(letters[1:length(xx[[r]])])
  }

  M <- list()

  for (t in 1:length(xx)) {
    aa <- xx[[t]]
    ll <- length(xx[[t]])

    M[[length(M) + 1]] <- ll

    # compute the average motif instance
    aa[["Mean"]] <- apply(aa[,1:ll], 1, mean)
    aa[["Mean"]] <- plyr::round_any(aa[["Mean"]], 0.01)

    # compute the quantile as shape
    aa[["Q80"]] <- apply(aa[,1:ll], 1, stats::quantile, probs =  0.8,
      na.rm = TRUE)
    aa[["Q80"]] <- plyr::round_any(aa[["Q80"]], 0.01)
    aa[["Median"]] <- apply(aa[,1:ll], 1, stats::quantile, probs =  0.5,
      na.rm = TRUE)
    aa[["Median"]] <- plyr::round_any(aa[["Median"]], 0.01)
    aa[["Q20"]] <- apply(aa[,1:ll], 1, stats::quantile, probs =  0.2,
      na.rm = TRUE)
    aa[["Q20"]] <- plyr::round_any(aa[["Q20"]], 0.01)

    # add the time as running index to the data frame
    aa[["Time"]] <- seq(1:nrow(aa))

    if (save == TRUE) {
      save(aa, file = paste0(path, top, "_aa_", t, "df.RData"))
    }

    normaleshape[[length(normaleshape) + 1]] <- aa
    names(normaleshape)[length(normaleshape)] <-
      paste0("Motif ", ifelse((t == length(xx) & no.noise == FALSE), "Noise", t))

    if (plot == TRUE) {
      # plot the result
      meltdf <- data.table::melt(aa[, c("Q80", "Q20", "Mean", "Median", "Time")],
        id = "Time")

      g <- ggplot2::ggplot(meltdf,
        ggplot2::aes(x = Time, y = value,
          colour = variable, group = variable)) +
        ggplot2::stat_summary(fun.data = "mean_cl_boot", geom = "smooth")

      grDevices::pdf(file = paste0(path, "Motif_",
        ifelse((t == length(xx) & no.noise == FALSE),
          "Noise", t), ".pdf"),
        height = 8,
        width = 12)

      print(g)

      grDevices::graphics.off()
    }
  }

  if (all == TRUE) {
    motifs <- lapply(normaleshape, function(x) x[[plot.var]])

    motifs <- as.data.frame(motifs)

    for (i in 1:length(names(motifs))) {
      names(motifs)[i] <- paste0(names(motifs)[i], ", n=", M[[i]])
    }

    motifs[["Time"]] <- seq(1, nrow(motifs))

    df <- data.table::melt(motifs,  id.vars = 'Time', variable.name = 'Load')

    p <- ggplot2::ggplot(df, ggplot2::aes(Time,value)) +
      ggplot2::geom_line(ggplot2::aes(colour = Load))

    q <- ggplot2::ggplot(df, ggplot2::aes(Time,value)) +
      ggplot2::geom_line() + ggplot2::facet_grid(Load ~ .)

    if (publication == TRUE) {
      p <- p + ggplot2::theme_light()
      p <- p + ggplot2::theme(legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 12),
        axis.text = ggplot2::element_text(size = 14),
        axis.title = ggplot2::element_text(size = 14,
          face = "bold"),
        legend.position = c(1,1),
        legend.background = ggplot2::element_rect(colour = "black"),
        legend.justification = c("right", "top"),
        legend.box.just = "top",
        legend.key = ggplot2::element_blank(),
        legend.text.align = 0,
        legend.margin = ggplot2::margin(t = 0.04, unit = "cm"),
        #legend.key.size = unit(0.45, "cm"),
        legend.key.height = ggplot2::unit(0.4,"cm"))

      q <- q + ggplot2::theme(legend.title = ggplot2::element_blank() ,
        legend.text = ggplot2::element_text(size = 14),
        axis.text = ggplot2::element_text(size = 14),
        axis.title = ggplot2::element_text(size = 14,
          face = "bold"))
    }

    pdf(file = paste0(path, "AllMotifs_", machine, ".pdf"),
      height = 8,
      width = 12,
      onefile = TRUE)

    print(p)
    print(q)

    graphics.off()
  }
  return(normaleshape)
}