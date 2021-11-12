# Script to extract profiles from motifs in Shahabs data set
# MIT License
# Copyright (c) 2021 Nicole Ludwig

library(data.table)
library(ggplot2)
library(Hmisc)

found.motifs = i = Time = value = variable = NULL

aggregation <- c("1 second", "2 min", "10 min", "2 hours")

for (a in aggregation) {
  load(paste0("data/motifs_", a, ".RData"))

  # get the number of motifs
  nr.motif <- length(found.motifs$Motif.raw)
  motif.names <- paste("Motif", 1:nr.motif)

  for (l in 1:nr.motif) {

    # get the number of instances in a motif
    nr.instances <- length(found.motifs$Motif.raw[[l]])

    # calculate the averages etc over all instances

    aa <- c()

    current <- as.data.frame(found.motifs$Motif.raw[[l]])
    names(current) <- letters[1:nr.instances]

    aa[["Mean"]] <- apply(current, 1, mean)
    aa[["Mean"]] <- plyr::round_any(aa[["Mean"]], 0.01)

    # compute the quantile as shape
    aa[["Q80"]] <- apply(current, 1, quantile, probs =  0.8,
                         na.rm = TRUE)
    aa[["Q80"]] <- plyr::round_any(aa[["Q80"]], 0.01)

    aa[["Median"]] <- apply(current, 1, quantile, probs =  0.5,
                            na.rm = TRUE)
    aa[["Median"]] <- plyr::round_any(aa[["Median"]], 0.01)

    aa[["Q20"]] <- apply(current, 1, quantile, probs =  0.2,
                         na.rm = TRUE)
    aa[["Q20"]] <- plyr::round_any(aa[["Q20"]], 0.01)

    # add the time as running index to the data frame
    aa[["Time"]] <- seq(1:length(aa[["Mean"]]))

    aa <- as.data.frame(aa)

    save(aa, file = paste0(a, "_aa_", l, "df.RData"))
    write.csv(aa, file = paste0(a, "_Motif_", l, "_data.csv"))

    # plot the result
    meltdf <- melt(aa[, c("Q80", "Q20", "Mean", "Median", "Time")], id = "Time")

    g <- ggplot(meltdf, aes(x = Time,
                            y = value,
                            colour = variable,
                            group = variable)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "smooth")

    pdf(file = paste0(a,"_Motif_", l, ".pdf"),
        height = 8,
        width = 12)

    print(g)
    graphics.off()
  }
}
