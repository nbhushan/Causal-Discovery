


#
# ###------- Load Results -------------------------------------------
# #first set source file
# files <- list.files(path = 'C:\\Users\\p275030\\iCloudDrive\\Documents\\Papers\\CD paper\\Final Analysis\\Results\\',pattern = '^cdSim[0-9]\\.rds$', full.names = TRUE)
#
# #BEST DAG (sid lower bound)
# res2 <- array(numeric(), c(2, 2, 4, 4, 4, 3, 3, 0))
# res2 <- getArray(readRDS("Results/cdSim00.rds"), "value")
#
# #Worst DAG (sid upper bound)
# res3 <- array(numeric(), c(2, 2, 4, 4, 4, 3, 3, 0))
# res3 <- getArray(readRDS("Results/cdSim01.rds"), "value")
#
#  FUNC <- function (x) getArray(readRDS(x), "value")
# res2 <- array(numeric(), c(2, 2, 4, 4, 4, 3, 3, 0))
# for (f in files) {
#   temp <- FUNC(f)
#   res2 <- abind(res2, temp, use.dnns = TRUE)
#   rm(temp)
# }
# rm(f, files)


correlateMeasures <- function(res2)
{
  names(dimnames(res2))[[1]] <- "metric"
  SHDRes <- res2[metric = "hamming", , , , , , ,]
  non.sim.margins <- setdiff(names(dimnames(SHDRes)), "n.sim")
  Var <- apply(SHDRes, non.sim.margins, mean)
  SHDwide <- spread(array2df(Var), key = algo, value = value)
  SIDRes <- res2[metric = "sid", , , , , , ,]
  non.sim.margins <- setdiff(names(dimnames(SIDRes)), "n.sim")
  Var <- apply(SIDRes, non.sim.margins, mean)
  SIDwide <- spread(array2df(Var), key = algo, value = value)
  PCmetric <- data.frame(hamming = SHDwide$pc, sid = SIDwide$pc)
  pcAgreement <- cor(x = rank(PCmetric$hamming),
                     y = rank(PCmetric$sid))
  lingamMetric <-
    data.frame(hamming = SHDwide$LINGAM, sid = SIDwide$LINGAM)
  lingamAgreement <-
    cor(x = rank(lingamMetric$hamming),
        y = rank(lingamMetric$sid))
  return(mean(pcAgreement, lingamAgreement))
}

getDf <- function(res2, metric = "hamming")
{
  names(dimnames(res2))[[1]] <- "metric"
  SIDRes <- res2[metric = metric, , , , , , ,]
  non.sim.margins <- setdiff(names(dimnames(SIDRes)), "n.sim")
  Var <- apply(SIDRes, non.sim.margins, mean)
  wide <- spread(array2df(Var), key = algo, value = value)
  return(wide)
}

se <- function(x)
  sqrt(var(x) / length(x))

getDf_se <- function(res2, metric = "hamming")
{
  names(dimnames(res2))[[1]] <- "metric"
  SIDRes <- res2[metric = metric, , , , , , ,]
  non.sim.margins <- setdiff(names(dimnames(SIDRes)), "n.sim")
  Var <- apply(SIDRes, non.sim.margins, se)
  wide <- spread(array2df(Var), key = algo, value = value)
  return(wide)
}

getDf_N <- function(res2, metric = "hamming")
{
  names(dimnames(res2))[[1]] <- "metric"
  resdf <- res2[metric = metric, , n = "300", , , , ,]
  non.sim.margins <- setdiff(names(dimnames(resdf)), "n.sim")
  Var <- apply(resdf, non.sim.margins, mean)
  wide <- spread(array2df(Var), key = algo, value = value)
  return(wide)
}

getDf_p <- function(res2, metric = "hamming")
{
  names(dimnames(res2))[[1]] <- "metric"
  resdf <- res2[metric = metric, , , p = "8", , , ,]
  non.sim.margins <- setdiff(names(dimnames(resdf)), "n.sim")
  Var <- apply(resdf, non.sim.margins, mean)
  wide <- spread(array2df(Var), key = algo, value = value)
  return(wide)
}


createGraph <-
  function(pd,
           metric = "hamming",
           xlabel = "",
           figName = "") {
    pdf(figName,
        paper = "a4r",
        width = 18,
        height = 15)
    
    if (metric == "hamming") {
      ylim <- c(-1.2, 1.2)
      ymin.refline = -1.1
      ymax.refline = -0.01
      ylab = "mean SHD"
      main = "Structural Hamming Distance"
    } else if (metric == "hammingSE") {
      ylim <- c(-0.02, 0.025)
      ymin.refline = -0.02
      ymax.refline = -0.001
      ylab = "se SHD"
      main = "Structural Hamming Distance"
    } else if (metric == "sid") {
      ylim <- c(-1.2, 1.2)
      ymin.refline = -1.1
      ymax.refline = -0.01
      ylab = "mean SID"
      main = "Structural Intervention Distance"
    } else {
      ylim <- c(-0.02, 0.025)
      ymin.refline = -0.02
      ymax.refline = -0.001
      ylab = "se SID"
      main = "Structural Intervention Distance"
    }
    
    
    plot(
      pd$LINGAM,
      #log=c("y"),
      type = "n",
      ylim = ylim,
      bty = "n",
      xlab = xlabel,
      ylab = ylab,
      las = 1,
      xaxt = "n"
    )
    
    #axis(side=2, at= c(0,1), labels=NULL)
    ##
    ## Add vertical lines (using R function lines.nestedloop)
    ##
    #lines(pd2, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE", "#FFFFF"))
    lines(pd, col = gray.colors(
      ncol(pd) - 2,
      start = 0.0,
      end = 0.9,
      alpha = 0.2
    ))
    ##
    ## Add reference lines (using R function lines.nestedloop)
    ##
    lines(
      pd,
      which = "r",
      ymin.refline = ymin.refline,
      ymax.refline = ymax.refline,
      cex.ref = 0.9
    )
    
    ##
    ## Estimates and legend (using standard R functions lines and legend)
    cols <- wes_palette("Darjeeling2", 2, "discrete")
    ##
    if (metric != "hamming" && metric != "hammingSE") {
      lines(pd$pc.lb, col = cols[1],  type = "s")
      lines(pd$pc.ub, col = cols[1],  type = "s") # pc
      polygon(
        x = c(1:nrow(pd), rev(1:nrow(pd))),
        y = c(pd$pc.lb, rev(pd$pc.ub)),
        border = NA,
        col = gray(0.9, alpha = 0.5)
      )
      lines(pd$LINGAM,
            col = cols[2],
            lwd = 1,
            type = "s")    # lingam
    } else{
      lines(pd$pc, col = cols[1],  type = "s")
      lines(pd$LINGAM,
            col = cols[2],
            lwd = 1,
            type = "s")
    }
    
    legend(
      "top",
      lwd = rep(2, 3),
      col = cols,
      cex = 0.8,
      c("PC", "LiNGAM"),
      horiz = T
    )
    ##
    dev.off()
  }
