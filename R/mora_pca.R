pca_fscores <- function(resPCA, design, axis1, axis2, col4obs, col4group, inference = TRUE){

   a.points <- 0.9
   if(inference){
      a.points <- 0.15
   }

   my.fi.plot <- PTCA4CATA::createFactorMap(resPCA$Fixed.Data$ExPosition.Data$fi,
                                 title = "Row Factor Scores",
                                 axis1 = axis1, axis2 = axis2,
                                 pch = 19,
                                 cex = 2,
                                 text.cex = 2.5,
                                 alpha.points = a.points,
                                 display.labels = FALSE,
                                 col.points = col4obs,
                                 col.labels = col4obs,
   )

   fi.labels <- PTCA4CATA::createxyLabels.gen(axis1, axis2,
                                   lambda = resPCA$Fixed.Data$ExPosition.Data$eigs,
                                   tau = round(resPCA$Fixed.Data$ExPosition.Data$t),
                                   axisName = "Component "
   )
   fi.plot <- my.fi.plot$zeMap + fi.labels

   if(inference){
      #update in a future version
      group.mean <- stats::aggregate(resPCA$Fixed.Data$ExPosition.Data$fi,
                              by = list(design), # must be a list
                              mean)

      # need to format the results from `aggregate` correctly
      rownames(group.mean) <- group.mean[,1] # Use the first column as row names
      fi.mean <- group.mean[,-1] # Exclude the first column

      fi.mean.plot <- PTCA4CATA::createFactorMap(fi.mean,
                                      alpha.points = 1,
                                      col.points = col4group[rownames(fi.mean)],
                                      col.labels = col4group[rownames(fi.mean)],
                                      pch = 17,
                                      cex = 3,
                                      text.cex = 4,
                                      axis1 = axis1,
                                      axis2 = axis2)
      fi.WithMean <- fi.plot + fi.mean.plot$zeMap_dots + fi.mean.plot$zeMap_text

      TI <- PTCA4CATA::MakeToleranceIntervals(resPCA$Fixed.Data$ExPosition.Data$fi,
                                   design = design,
                                   col = col4group[rownames(fi.mean)],
                                   p.level = 0.95,
                                   axis1 = axis1,
                                   axis2 = axis2)
      fi.WithTI <- fi.WithMean + TI
      print(fi.WithTI)


      fi.boot <- PTCA4CATA::Boot4Mean(resPCA$Fixed.Data$ExPosition.Data$fi,
                                       design = design,
                                       niter = 1000, suppressProgressBar = FALSE)


      bootCI4mean <- PTCA4CATA::MakeCIEllipses(fi.boot$BootCube,
                                    col = col4group[rownames(fi.mean)],
                                    p.level = 0.95,
                                    axis1 = axis1,
                                    axis2 = axis2,
                                    )

      fi.WithMeanCI <- fi.WithMean + bootCI4mean
      print(fi.WithMeanCI)
   }
   else{
      print(fi.plot)
   }

}

#' @import ggplot2
pca_columns <- function(resPCA, data, axis1, axis2, col4var = NULL, important = FALSE){

   cor.loading <- stats::cor(data, resPCA$Fixed.Data$ExPosition.Data$fi)

   colFactorMap <- col4var
   if(is.null(col4var)){
      colFactorMap <- "darkorchid4"
   }

   loading.plot <- PTCA4CATA::createFactorMap(cor.loading,
                                   constraints = list(minx = -1, miny = -1, maxx = 1, maxy = 1),
                                   font.face = "plain",
                                   text.cex = 2.5,
                                   cex = 1,
                                   force = 0.5,
                                   axis1 = axis1,
                                   axis2 = axis2,
                                   col.points = colFactorMap,
                                   col.labels = colFactorMap)

   LoadingMapWithCircles <- loading.plot$zeMap +
      PTCA4CATA::addArrows(cor.loading,
                color = "black",
                size = 0.5,
                arrowLength = 0.2,
                alpha = 0.2,
                axis1 = axis1,
                axis2 = axis2) +
      PTCA4CATA::addCircleOfCor() + xlab(paste("Component", axis1)) + ylab(paste("Component", axis2))

   print(LoadingMapWithCircles)



   my.fj.plot <- PTCA4CATA::createFactorMap(resPCA$Fixed.Data$ExPosition.Data$fj,
                              font.face = "plain",
                              text.cex = 2.5,
                              cex = 1,
                              force = 0.5,
                              axis1 = axis1,
                              axis2 = axis2,
                              col.points = colFactorMap,
                              col.labels = colFactorMap
   )

   fj.labels <- PTCA4CATA::createxyLabels.gen(axis1, axis2,
                                   lambda = resPCA$Fixed.Data$ExPosition.Data$eigs,
                                   tau = round(resPCA$Fixed.Data$ExPosition.Data$t),
                                   axisName = "Component "
   )
   fj.plot <- my.fj.plot$zeMap + fj.labels
   print(fj.plot)

   signed.ctrJ <- resPCA$Fixed.Data$ExPosition.Data$cj *
      sign(resPCA$Fixed.Data$ExPosition.Data$fj)

   # plot contributions for component 1
   ctrJ.1 <- PTCA4CATA::PrettyBarPlot2(signed.ctrJ[,axis1],
                            threshold = 1 / NROW(signed.ctrJ),
                            font.size = 3,
                            color4bar = col4var,
                            ylab = 'Contributions',
                            ylim = c(1.2*min(signed.ctrJ[,axis1:axis2]), 1.2*max(signed.ctrJ[,axis1:axis2])),
                            signifOnly = important) +
      ggtitle("Contribution barplots", subtitle = paste("Component", axis1))

   # plot contributions for component 2
   ctrJ.2 <- PTCA4CATA::PrettyBarPlot2(signed.ctrJ[,axis2],
                            threshold = 1 / NROW(signed.ctrJ),
                            font.size = 3,
                            color4bar = col4var,
                            ylab = 'Contributions',
                            ylim = c(1.2*min(signed.ctrJ[,axis1:axis2]), 1.2*max(signed.ctrJ[,axis1:axis2])),
                            signifOnly = important) +
      ggtitle("",subtitle = paste("Component", axis2))

   print(ctrJ.1)
   print(ctrJ.2)

   BR <- resPCA$Inference.Data$fj.boots$tests$boot.ratios

   ba001.BR1 <- PTCA4CATA::PrettyBarPlot2(BR[,axis1],
                               threshold = 2,
                               font.size = 3,
                               color4bar = col4var,
                               ylab = 'Bootstrap ratios',
                               ylim = c(1.2*min(BR[,axis1:axis2]), 1.2*max(BR[,axis1:axis2])),
                               signifOnly = important)+
      ggtitle("Bootstrap ratios", subtitle = paste0('Component ', axis1))

   ba002.BR2 <- PTCA4CATA::PrettyBarPlot2(BR[,axis2],
                                          threshold = 2,
                                          font.size = 3,
                                          color4bar = col4var,
                                          ylab = 'Bootstrap ratios',
                                          ylim = c(1.2*min(BR[,axis1:axis2]), 1.2*max(BR[,axis1:axis2])),
                                          signifOnly = important)+
      ggtitle("Bootstrap ratios", subtitle = paste0('Component ', axis2))

   print(ba001.BR1)
   print(ba002.BR2)

   pca_columns_res <- list(corCircle = LoadingMapWithCircles,
                           fj.plot = my.fj.plot,
                           ctr = list(ctrJ.1, ctrJ.2),
                           BR = list(ba001.BR1, ba002.BR2))

   invisible(pca_columns_res)

}

#' PCA with inference and graphs.
#' @param data The numerical data for PCA.
#' @param design (default = NULL) A vector or dummy-coded matrix that describes the rows.
#' @param make_design_nominal (default = TRUE) If design is a dummy-coded matrix, should be set to FALSE.
#' @param col4obs (default = "olivedrab3") A single color or vector of colors whose length is equal to nrow(data).
#' @param col4group (default = "olivedrab3") A single color or vector of colors whose length is the number of groups in design.
#' @param col4var (default = NULL) A vector of colors whose length is the equal to ncol(data).
#' @param center (default = TRUE) Whether to center variables
#' @param scale (default = "SS1") Whether to scale variables
#' @param want34 (default = FALSE) By default, only prints dimensions 1 and 2. Set to TRUE for 3 and 4.
#' @param inference (default = TRUE) When design contains 2 or more groups, computes CI and TI on observations. FALSE if no groups.
#' @param important (default = FALSE) Only display variables above threshold for ctr and BR
#' @param test_iters (default = 1000) How many permutations to run
#' @export
mora_pca <- function(data,
                    design = NULL,
                    make_design_nominal = TRUE,
                    col4obs = "olivedrab3",
                    col4group = "olivedrab3",
                    col4var = NULL,
                    center = TRUE,
                    scale = "SS1",
                    want34 = FALSE,
                    inference = TRUE,
                    important = FALSE,
                    test_iters = 1000){
   data_cor <- stats::cor(data)

   corrplot::corrplot(data_cor, tl.cex = 0.7, tl.pos = "lt", tl.col = "black",
            addCoefasPercent = TRUE, addCoef.col = "black",
            number.cex = 0.5, method = "color")

   resPCA <- InPosition::epPCA.inference.battery(data, center = center, scale = scale,
                                     DESIGN = design, graphs = FALSE, test.iters = test_iters,
                                     make_design_nominal = make_design_nominal)

   scree <- PTCA4CATA::PlotScree(ev = resPCA$Fixed.Data$ExPosition.Data$eigs,
                           p.ev = resPCA$Inference.Data$components$p.vals,
                           plotKaiser = TRUE)

   pca_fscores(resPCA = resPCA, design = design, axis1 = 1, axis2 = 2,
               col4obs = col4obs, col4group = col4group, inference = inference)

   pca_columns(resPCA, data = data, 1, 2, col4var = col4var, important = important)

   if(want34){
      pca_fscores(resPCA = resPCA, design = design, axis1 = 3, axis2 = 4,
                  col4obs = col4obs, col4group = col4group, inference = inference)

      pca_columns(resPCA, data = data, 3, 4, col4var = col4var, important = important)
   }
}
