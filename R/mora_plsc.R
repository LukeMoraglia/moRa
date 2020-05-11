#' @import prettyGraphs
#' @import ggplot2
#' @importFrom ggplotify as.grob
#' @importFrom gridExtra grid.arrange
plsc_saliences <- function(resPLS, lvNum = 1, important = FALSE){
   p <- resPLS$TExPosition.Data$pdq$p
   q <- resPLS$TExPosition.Data$pdq$q



   plot.p <- PrettyBarPlot2(p[,lvNum],
                            threshold = sqrt(1/nrow(p)),
                            font.size = 3,
                            ylab = 'p for Lx',
                            horizontal = FALSE,
                            signifOnly = important
   ) +
      ggtitle(paste("Column saliences of X for Latent Variable", lvNum))

   plot.q <- PrettyBarPlot2(q[,lvNum],
                            threshold = sqrt(1/nrow(q)),
                            font.size = 3,
                            ylab = 'q for Ly',
                            signifOnly = important
   ) +
      ggtitle(paste("Column saliences of Y for Latent Variable", lvNum))

   pPlot <- length(plot.p$data$bootratio) > 0
   qPlot <- length(plot.q$data$bootratio) > 0

   if(pPlot & qPlot){
       grid.arrange(as.grob(plot.p), as.grob(plot.q), nrow = 2)
   }
   else if(pPlot){
      print(plot.p)
   }
   else if(qPlot){
      print(plot.q)
   }

}

#' @import data4PCCAR
#' @import prettyGraphs
#' @import ggplot2
#' @importFrom ggplotify as.grob
#' @importFrom gridExtra grid.arrange
plsc_boot_ratio <- function(data1, data2,
                    center1, center2,
                    scale1, scale2,
                    Fi = NULL, Fj = NULL, lvNum = 1,
                    important = FALSE){


      boot.res <- Boot4PLSC(data1, data2, center1, center2, scale1, scale2, Fi, Fj)

      BRIJ <- firstpos(boot.res$bootRatios.i, boot.res$bootRatios.j)
      boot.res$bootRatios.i <- BRIJ$P
      boot.res$bootRatios.j <- BRIJ$Q

      plotBRi <- PrettyBarPlot2(boot.res$bootRatios.i[,lvNum],
                                threshold = 2,
                                font.size = 3,
                                signifOnly = important,
                                ylab = 'BR for Lx',
                                horizontal = FALSE
                                ) +
         ggtitle(paste("Bootstrap ratios of X for Latent Variable", lvNum))

      plotBRj <- PrettyBarPlot2(boot.res$bootRatios.j[,lvNum],
                                threshold = 2,
                                font.size = 3,
                                signifOnly = important,
                                ylab = 'BR for Ly',
                                horizontal = TRUE
      ) +
         ggtitle(paste("Bootstrap ratios of Y for Latent Variable", lvNum))

      briplot <- length(plotBRi$data$bootratio) > 0
      brjplot <- length(plotBRj$data$bootratio) > 0

      if(briplot & brjplot){
          grid.arrange(as.grob(plotBRi), as.grob(plotBRj), nrow = 2)
      }
      else if(briplot){
         print(plotBRi)
      }
      else if(brjplot){
         print(plotBRj)
      }
}

#' @import PTCA4CATA
#' @import stats
#' @import ggplot2
#' @import ExPosition
plsc_latent_variable <- function(resPLS, design, col4obs, col4group, lvNum, inference = inference){
   latvar <- cbind(resPLS$TExPosition.Data$lx[,lvNum],resPLS$TExPosition.Data$ly[,lvNum])
   colnames(latvar) <- c(paste("Lx", lvNum), paste("Ly", lvNum))

   a.points <- 0.9
   if(inference){
      a.points <- 0.2
   }

   plot.lv <- createFactorMap(latvar,
                              col.points = col4obs,
                              col.labels = col4obs,
                              alpha.points = a.points,
                              display.labels = FALSE
   )

   lat.cor <- cor(latvar[,1], latvar[,2])

   lv.label <- labs(x = paste0("Singular value = ", round(resPLS$TExPosition.Data$pdq$Dv[lvNum], 3), " Correlation = ",
                               round(lat.cor, 3)))

   if(inference){
      # compute means
      lv.group <- getMeans(latvar, design)

      # get bootstrap intervals of groups
      lv.group.boot <- Boot4Mean(latvar, design, niter = 1000)
      colnames(lv.group.boot$BootCube) <- c(paste("Lx", lvNum), paste("Ly", lvNum))


      plot.mean <- createFactorMap(lv.group,
                                    col.points = col4group[rownames(lv.group)],
                                    col.labels = col4group[rownames(lv.group)],
                                    cex = 4,
                                    pch = 17,
                                    alpha.points = 0.8)

      plot.meanCI <- MakeCIEllipses(lv.group.boot$BootCube[,c(1:2),],
                                     col = col4group[rownames(lv.group)],
                                     names.of.factors = c(paste("Lx", lvNum), paste("Ly", lvNum))
      )

      plot<- plot.lv$zeMap_background + plot.lv$zeMap_dots +
         plot.mean$zeMap_dots + plot.mean$zeMap_text + plot.meanCI + lv.label
      print(plot)
   }
   else{
      fullPlot <- plot.lv$zeMap + lv.label
      print(fullPlot)
   }
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
#' @import stats
#' @import corrplot
#' @import TExPosition
#' @import PTCA4CATA
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
                    important = FALSE){

   data_cor <- cor(data1, data2)

   corrplot(data_cor, tl.cex = 0.7, tl.pos = "lt", tl.col = "black",
            addCoefasPercent = TRUE, addCoef.col = "black",
            number.cex = 0.5, method = "color")

   resPLS <- tepPLS(data1, data2, scale1 = scale1, scale2 = scale2,
                    center1 = center1, center2 = center2,
                    DESIGN = design, graphs = FALSE)

   resPLS <- plsc_first_pos(resPLS, data1, data2, center1, center2, scale1, scale2)

   p <- resPLS$TExPosition.Data$pdq$p
   Dv <- resPLS$TExPosition.Data$pdq$Dv
   q <- resPLS$TExPosition.Data$pdq$q

   PLSperm <- perm4PLSC(data1, data2, scale1 = scale1, scale2 = scale2,
                        center1 = center1, center2 = center2,
                        permType = 'byColumns')

   singularScree <- PlotScree(resPLS$TExPosition.Data$pdq$Dv,
                              title = "Singular Values")

   eigenScree <- PlotScree(resPLS$TExPosition.Data$eigs,
                           p.ev = PLSperm$pEigenvalues,
                           plotKaiser = TRUE,
                           title = "Eigenvalues")

   plsc_latent_variable(resPLS, design, col4obs, col4group, 1, inference = inference)

   plsc_saliences(resPLS, 1, important)
   plsc_boot_ratio(data1, data2, center1, center2, scale1, scale2, resPLS$TExPosition.Data$fi,
           resPLS$TExPosition.Data$fj, 1, important)

   deflation1 <- (p[,1] * Dv[1]) %*% t(q[,1])
   data_cor2 <- data_cor - deflation1

   corrplot(data_cor2, tl.cex = 0.7, tl.pos = "lt", tl.col = "black",
            addCoefasPercent = TRUE, addCoef.col = "black",
            number.cex = 0.5, method = "color")

   plsc_latent_variable(resPLS, design, col4obs, col4group, 2, inference = inference)
   plsc_saliences(resPLS, 2, important)
   plsc_boot_ratio(data1, data2, center1, center2, scale1, scale2, resPLS$TExPosition.Data$fi,
           resPLS$TExPosition.Data$fj, 2, important)

}

#' @import data4PCCAR
plsc_first_pos <- function(resFile, data1, data2, center1, center2, scale1, scale2){
   if (attr(resFile, "class")[1] != "texpoOutput") {
      stop("plsc_first_pos works only with TExPosition objects")
   }

   data1.CS <- expo.scale(data1, center1, scale1)
   data2.CS <- expo.scale(data2, center2, scale2)

   QP <- firstpos(resFile$TExPosition.Data$pdq$p, resFile$TExPosition.Data$pdq$q)
   resFile$TExPosition.Data$pdq$p <- QP$P
   resFile$TExPosition.Data$pdq$q <- QP$Q

   resFile$TExPosition.Data$lx <- data1.CS %*% resFile$TExPosition.Data$pdq$p
   resFile$TExPosition.Data$ly <- data2.CS %*% resFile$TExPosition.Data$pdq$q

   return(resFile)
}
