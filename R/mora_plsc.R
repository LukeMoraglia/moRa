#' @import ggplot2
plsc_saliences <- function(resPLS, lvNum = 1, important = FALSE){
   p <- resPLS$TExPosition.Data$pdq$p
   q <- resPLS$TExPosition.Data$pdq$q



   plot.p <- PTCA4CATA::PrettyBarPlot2(p[,lvNum],
                            threshold = sqrt(1/nrow(p)),
                            font.size = 3,
                            ylim = c(min(-sqrt(1/nrow(p)), p[,lvNum]),
                                     max(sqrt(1/nrow(p)), p[,lvNum])),
                            ylab = 'p for Lx',
                            horizontal = FALSE,
                            signifOnly = important
   ) +
      ggtitle(paste("Column saliences of X for Latent Variable", lvNum))

   plot.q <- PTCA4CATA::PrettyBarPlot2(q[,lvNum],
                            threshold = sqrt(1/nrow(q)),
                            font.size = 3,
                            ylab = 'q for Ly',
                            signifOnly = important
   ) +
      ggtitle(paste("Column saliences of Y for Latent Variable", lvNum))

   pPlot <- length(plot.p$data$bootratio) > 0
   qPlot <- length(plot.q$data$bootratio) > 0

   results <- list()

   if(pPlot & qPlot){
       gridExtra::grid.arrange(ggplotify::as.grob(plot.p),
                               ggplotify::as.grob(plot.q),
                               nrow = 2)
      results <- list(pPlot = pPlot,
                      qPlot = qPlot)
   }
   else if(pPlot){
      print(plot.p)
      results <- list(pPlot = pPlot)
   }
   else if(qPlot){
      print(plot.q)
      results <- list(qPlot = qPlot)
   }

   invisible(results)

}


#' @import ggplot2
plsc_boot_ratio <- function(data1, data2,
                    center1, center2,
                    scale1, scale2,
                    Fi = NULL, Fj = NULL, lvNum = 1,
                    important = FALSE){


      boot.res <- data4PCCAR::Boot4PLSC(data1, data2, center1, center2,
                                        scale1, scale2, Fi, Fj)

      BRIJ <- data4PCCAR::firstpos(boot.res$bootRatios.i, boot.res$bootRatios.j)
      boot.res$bootRatios.i <- BRIJ$P
      boot.res$bootRatios.j <- BRIJ$Q

      plotBRi <- PTCA4CATA::PrettyBarPlot2(boot.res$bootRatios.i[,lvNum],
                                threshold = 2,
                                font.size = 3,
                                ylim = c(min(-2, 1.2*min(boot.res$bootRatios.i[,lvNum])),
                                         max(2, 1.2*max(boot.res$bootRatios.i[,lvNum]))),
                                signifOnly = important,
                                ylab = 'BR for Lx',
                                horizontal = FALSE
                                ) +
         ggtitle(paste("Bootstrap ratios of X for Latent Variable", lvNum))

      plotBRj <- PTCA4CATA::PrettyBarPlot2(boot.res$bootRatios.j[,lvNum],
                                threshold = 2,
                                font.size = 3,
                                ylim = c(min(-2, 1.2*min(boot.res$bootRatios.j[,lvNum])),
                                         max(2, 1.2*max(boot.res$bootRatios.j[,lvNum]))),
                                signifOnly = important,
                                ylab = 'BR for Ly',
                                horizontal = TRUE
      ) +
         ggtitle(paste("Bootstrap ratios of Y for Latent Variable", lvNum))

      briplot <- length(plotBRi$data$bootratio) > 0
      brjplot <- length(plotBRj$data$bootratio) > 0

      results <- list()
      if(briplot & brjplot){
         gridExtra::grid.arrange(ggplotify::as.grob(plotBRi),
                                 ggplotify::as.grob(plotBRj),
                                 nrow = 2)
         results <- list(BRx = plotBRi,
                         BRy = plotBRj)
      }
      else if(briplot){
         print(plotBRi)
         results <- list(BRx = plotBRi)
      }
      else if(brjplot){
         print(plotBRj)
         results <- list(BRy = plotBRj)
      }

      invisible(results)
}


#' @import ggplot2
plsc_latent_variable <- function(resPLS,
                                 design,
                                 col4obs,
                                 col4group,
                                 lvNum,
                                 inference = TRUE){

   latvar <- cbind(resPLS$TExPosition.Data$lx[,lvNum],resPLS$TExPosition.Data$ly[,lvNum])
   colnames(latvar) <- c(paste("Lx", lvNum), paste("Ly", lvNum))

   a.points <- 0.9
   if(inference){
      a.points <- 0.2
   }

   plot.lv <- PTCA4CATA::createFactorMap(latvar,
                              col.points = col4obs,
                              col.labels = col4obs,
                              alpha.points = a.points,
                              display.labels = FALSE
   )

   lat.cor <- stats::cor(latvar[,1], latvar[,2])

   lv.label <- labs(x = paste0("Singular value = ",
                               round(resPLS$TExPosition.Data$pdq$Dv[lvNum], 3),
                               " Correlation = ",
                               round(lat.cor, 3)))

   if(inference){
      # compute means
      lv.group <- PTCA4CATA::getMeans(latvar, design)

      # get bootstrap intervals of groups
      lv.group.boot <- PTCA4CATA::Boot4Mean(latvar, design, niter = 1000)
      colnames(lv.group.boot$BootCube) <- c(paste("Lx", lvNum), paste("Ly", lvNum))


      plot.mean <- PTCA4CATA::createFactorMap(lv.group,
                                    col.points = col4group[rownames(lv.group)],
                                    col.labels = col4group[rownames(lv.group)],
                                    cex = 4,
                                    pch = 17,
                                    alpha.points = 0.8)

      plot.meanCI <- PTCA4CATA::MakeCIEllipses(lv.group.boot$BootCube[,c(1:2),],
                                     col = col4group[rownames(lv.group)],
                                     names.of.factors = c(paste("Lx", lvNum),
                                                          paste("Ly", lvNum))
      )

      plot<- plot.lv$zeMap_background + plot.lv$zeMap_dots +
         plot.mean$zeMap_dots + plot.mean$zeMap_text + plot.meanCI + lv.label
      print(plot)
   }
   else{
      plot <- plot.lv$zeMap + lv.label
      print(plot)
   }

   results <- list(lvplot = plot)
   invisible(results)
}

#' PLSC with inferences and graphs.
#' @param data1 A numerical data table with observations on rows and variables on columns.
#' @param data2 A numerical data table with the same obs on row and different vars on columns.
#' @param design (default = NULL) A design vector or matrix for the rows
#' @param make_design_nominal (default = TRUE) If TRUE, design is a vector. If FALSE, design is a dummy-coded matrix.
#' @param col4obs (default = "olivedrab3") A single color or vector of colors whose length is equal to nrow(data1).
#' @param col4group (default = "olivedrab3") A single color or vector of colors whose length is the number of groups in design.
#' @param center1 (default = TRUE) Whether to center variables in data1.
#' @param center2 (default = TRUE) Whether to center variables in data2.
#' @param scale1 (default = "SS1") Whether to scale variables in data1.
#' @param scale2 (default = "SS1") Whether to scale variables in data2.
#' @param inference (default = TRUE) When design contains 2 or more groups, computes CI and TI on observations. FALSE if no groups.
#' @param important (default = FALSE) If TRUE, graphs have only the important saliences/bootstrap ratios.
#' @param corrplot (default = TRUE) Whether to show a correlation matrix plot
#' @export
mora_plsc <- function(data1,
                    data2,
                    design = NULL,
                    make_design_nominal = TRUE,
                    col4obs = "olivedrab3",
                    col4group = "olivedrab3",
                    center1 = TRUE,
                    center2 = TRUE,
                    scale1 = "SS1",
                    scale2 = "SS1",
                    inference = TRUE,
                    important = FALSE,
                    corrplot = TRUE){

   data_cor <- stats::cor(data1, data2)

   if(corrplot){
      corrplot::corrplot(data_cor, tl.cex = 0.7, tl.pos = "lt", tl.col = "black",
               addCoefasPercent = TRUE, addCoef.col = "black",
               number.cex = 0.5, method = "color")
   }

   resPLS <- TExPosition::tepPLS(data1, data2, scale1 = scale1, scale2 = scale2,
                    center1 = center1, center2 = center2,
                    DESIGN = design, graphs = FALSE)

   resPLS <- plsc_first_pos(resPLS, data1, data2, center1, center2, scale1, scale2)

   p <- resPLS$TExPosition.Data$pdq$p
   Dv <- resPLS$TExPosition.Data$pdq$Dv
   q <- resPLS$TExPosition.Data$pdq$q

   PLSperm <- data4PCCAR::perm4PLSC(data1, data2, scale1 = scale1, scale2 = scale2,
                        center1 = center1, center2 = center2,
                        permType = 'byColumns')

   # singularScree <- PTCA4CATA::PlotScree(resPLS$TExPosition.Data$pdq$Dv,
   #                            title = "Singular Values")

   eigenScree <- PTCA4CATA::PlotScree(resPLS$TExPosition.Data$eigs,
                           p.ev = PLSperm$pEigenvalues,
                           plotKaiser = TRUE,
                           title = "Eigenvalues")
   scree <- grDevices::recordPlot()

   lv1 <- plsc_latent_variable(resPLS, design, col4obs, col4group, 1, inference = inference)

   sal1 <- plsc_saliences(resPLS, 1, important)
   br1 <- plsc_boot_ratio(data1, data2, center1, center2, scale1, scale2, resPLS$TExPosition.Data$fi,
           resPLS$TExPosition.Data$fj, 1, important)

   if(corrplot){
      deflation1 <- (p[,1] * Dv[1]) %*% t(q[,1])
      data_cor2 <- data_cor - deflation1

      corrplot::corrplot(data_cor2, tl.cex = 0.7, tl.pos = "lt", tl.col = "black",
               addCoefasPercent = TRUE, addCoef.col = "black",
               number.cex = 0.5, method = "color")
   }
   lv2 <- plsc_latent_variable(resPLS, design, col4obs, col4group, 2, inference = inference)
   sal2 <- plsc_saliences(resPLS, 2, important)
   br2 <- plsc_boot_ratio(data1, data2, center1, center2, scale1, scale2, resPLS$TExPosition.Data$fi,
           resPLS$TExPosition.Data$fj, 2, important)

   results <- list(TExPosition.Data = resPLS$TExPosition.Data,
                   Plotting.Data = resPLS$Plotting.Data,
                   perm = PLSperm,
                   scree = scree,
                   lv1 = lv1,
                   lv2 = lv2,
                   sal1 = sal1,
                   sal2 = sal2,
                   br1 = br1,
                   br2 = br2
                   )
   invisible(results)
}


plsc_first_pos <- function(resFile, data1, data2, center1, center2, scale1, scale2){
   # if (attr(resFile, "class")[1] != "texpoOutput") {
   #    stop("plsc_first_pos works only with TExPosition objects")
   # }

   data1.CS <- ExPosition::expo.scale(data1, center1, scale1)
   data2.CS <- ExPosition::expo.scale(data2, center2, scale2)

   QP <- data4PCCAR::firstpos(resFile$TExPosition.Data$pdq$p, resFile$TExPosition.Data$pdq$q)
   resFile$TExPosition.Data$pdq$p <- QP$P
   resFile$TExPosition.Data$pdq$q <- QP$Q

   resFile$TExPosition.Data$lx <- data1.CS %*% resFile$TExPosition.Data$pdq$p
   resFile$TExPosition.Data$ly <- data2.CS %*% resFile$TExPosition.Data$pdq$q

   return(resFile)
}
