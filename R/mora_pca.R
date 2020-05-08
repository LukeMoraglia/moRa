#' @import ggplot2
#' @import stats
MakeToleranceIntervals1 <- function (data, design, axis1 = 1, axis2 = 2, names.of.factors = paste0("Dimension ",
                                                                                                   c(axis1, axis2)), col = NULL, centers = NULL, line.size = 1,
                                     line.type = 1, alpha.ellipse = 0.3, alpha.line = 0.5, p.level = 0.66,
                                     type = "hull")
{
   design <- factor(design)
   Nom2Rows <- levels(design)
   X <- data[, c(axis1, axis2)]
   if (length(design) != NROW(X)) {
      stop("Length of Design should be equal to nrow(Data)")
   }
   if (is.null(names.of.factors)) {
      names.of.factors = unlist(dimnames(X)[2])
   }
   if (is.null(names.of.factors)) {
      names.of.factors = paste0("Dimension ", c(axis1,
                                                axis2))
   }
   colnames(X) <- names.of.factors
   nItems = length(Nom2Rows)
   if (is.null(col)) {
      items.colors <- prettyGraphs::prettyGraphsColorSelection(nItems)
   }
   else {
      items.colors <- col
   }
   if (length(items.colors) == 1) {
      items.colors = rep(items.colors, nItems)
   }
   if (length(items.colors) != nItems) {
      items.colors = rep(items.colors[1], nItems)
   }
   LeGraph.elli <- list()
   for (i in 1:nItems) {
      X2plot <- as.data.frame(X[design == Nom2Rows[i], ])
      if (!is.null(centers)) {
         truc <- as.data.frame(sweep(sweep(as.matrix(X2plot),
                                           2, colMeans(X2plot)), 2, -as.matrix(centers[i,
                                                                                       ])))
         X2plot <- truc
      }
      if (tolower(type) == "ellipse") {
         elli <- ggplot2::stat_ellipse(data = X2plot[, c(axis1,
                                                         axis2)], ggplot2::aes_(color = alpha(items.colors[i],
                                                                                              alpha.line)), show.legend = FALSE, geom = "polygon",
                                       fill = ggplot2::alpha(items.colors[i], alpha.ellipse),
                                       type = "t", na.rm = TRUE, level = p.level,
                                       color = items.colors[i], size = line.size, linetype = line.type)
      }
      else {
         row2keep <- stats::complete.cases(X2plot[, 1:2])
         X.non.na <- X2plot[row2keep, 1:2]
         elli <- ggConvexHull(data = X.non.na, x_axis = 1,
                              y_axis = 2, percentage = p.level, col.line = items.colors[i],
                              alpha.line = alpha.line, line.size = line.size,
                              line.type = line.type, col.hull = items.colors[i],
                              alpha.hull = alpha.ellipse, names.of.factors = names.of.factors)
      }
      LeGraph.elli[[i]] <- elli
   }
   return(LeGraph.elli)
}


#' @import ggplot2
#' @import prettyGraphs
MakeCIEllipses1 <- function (data, axis1 = 1, axis2 = 2, names.of.factors = paste0("Dimension ",
                                                                                   c(axis1, axis2)), col = NULL, centers = NULL, line.size = 1,
                             line.type = 1, alpha.ellipse = 0.3, alpha.line = 0.5, p.level = 0.95)
{
   Nom2Rows <- unlist(dimnames(data)[1])
   X <- aperm(data, c(1, 3, 2))
   if (is.null(names.of.factors)) {
      names.of.factors = unlist(dimnames(data)[2])
   }
   DimBoot <- dim(data)
   dim(X) <- c(DimBoot[1] * DimBoot[3], DimBoot[2])
   rownames(X) <- rep(Nom2Rows, DimBoot[3])
   colnames(X) <- names.of.factors
   nItems = DimBoot[1]
   if (is.null(col)) {
      items.colors <- prettyGraphs::prettyGraphsColorSelection(nItems)
   }
   else {
      items.colors <- col
   }
   if (length(items.colors) == 1) {
      items.colors = rep(items.colors, nItems)
   }
   if (length(items.colors) != nItems) {
      items.colors = rep(items.colors[1], nItems)
   }
   LeGraph.elli <- list()
   for (i in 1:nItems) {
      X2plot <- as.data.frame(X[row.names(X) == Nom2Rows[i],
                                ])
      if (!is.null(centers)) {
         truc <- as.data.frame(sweep(sweep(as.matrix(X2plot),
                                           2, colMeans(X2plot)), 2, -as.matrix(centers[i,
                                                                                       ])))
         X2plot <- truc
      }
      elli <- ggplot2::stat_ellipse(data = X2plot[, 1:2], ggplot2::aes(color = alpha(items.colors[i],
                                                                                          alpha.line)), show.legend = FALSE, geom = "polygon",
                                    fill = ggplot2::alpha(items.colors[i], alpha.ellipse),
                                    type = "t", level = p.level, color = items.colors[i],
                                    size = line.size, linetype = line.type)
      LeGraph.elli[[i]] <- elli
   }
   return(LeGraph.elli)
}

#' @import PTCA4CATA
#' @import stats
pca_fscores <- function(resPCA, design, axis1, axis2, col4obs, col4group, inference = TRUE){
   #resPCA$Fixed.Data$Plotting.Data$fi.col <- col4fii
   a.points <- 0.9
   if(inference){
      a.points <- 0.15
   }

   my.fi.plot <- createFactorMap(resPCA$Fixed.Data$ExPosition.Data$fi,
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

   fi.labels <- createxyLabels.gen(axis1, axis2,
                                   lambda = resPCA$Fixed.Data$ExPosition.Data$eigs,
                                   tau = round(resPCA$Fixed.Data$ExPosition.Data$t),
                                   axisName = "Component "
   )
   fi.plot <- my.fi.plot$zeMap + fi.labels

   if(inference){
      group.mean <- aggregate(resPCA$Fixed.Data$ExPosition.Data$fi,
                              by = list(design), # must be a list
                              mean)

      # need to format the results from `aggregate` correctly
      rownames(group.mean) <- group.mean[,1] # Use the first column as row names
      fi.mean <- group.mean[,-1] # Exclude the first column

      fi.mean.plot <- createFactorMap(fi.mean,
                                      alpha.points = 1,
                                      col.points = col4group[rownames(fi.mean)],
                                      col.labels = col4group[rownames(fi.mean)],
                                      pch = 17,
                                      cex = 3,
                                      text.cex = 4,
                                      axis1 = axis1,
                                      axis2 = axis2)
      fi.WithMean <- fi.plot + fi.mean.plot$zeMap_dots + fi.mean.plot$zeMap_text

      TI <- MakeToleranceIntervals1(resPCA$Fixed.Data$ExPosition.Data$fi,
                                   design = design,
                                   col = col4group[rownames(fi.mean)],
                                   p.level = 0.95,
                                   axis1 = axis1,
                                   axis2 = axis2)
      fi.WithTI <- fi.WithMean + TI
      print(fi.WithTI)

      # Depending on the size of your data, this might take a while
      fi.boot <- Boot4Mean(resPCA$Fixed.Data$ExPosition.Data$fi,
                           design = design,
                           niter = 1000, suppressProgressBar = FALSE)


      bootCI4mean <- MakeCIEllipses1(fi.boot$BootCube[,c(axis1, axis2),],
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

#' @import stats
#' @import PTCA4CATA
#' @import ggplot2
#' @import prettyGraphs
pca_columns <- function(resPCA, data, axis1, axis2){
   #resPCA$Fixed.Data$Plotting.Data$fj.col <- col4Xvar
   cor.loading <- cor(data, resPCA$Fixed.Data$ExPosition.Data$fi)

   loading.plot <- createFactorMap(cor.loading,
                                   constraints = list(minx = -1, miny = -1, maxx = 1, maxy = 1),
                                   font.face = "plain",
                                   text.cex = 2.5,
                                   cex = 1,
                                   force = 0.5,
                                   axis1 = axis1,
                                   axis2 = axis2
                                   #col.points = res_pcaInf$Fixed.Data$Plotting.Data$fj.col,
                                   #col.labels = res_pcaInf$Fixed.Data$Plotting.Data$fj.col
   )

   LoadingMapWithCircles <- loading.plot$zeMap +
      addArrows(cor.loading,
                color = "black",
                size = 0.5,
                arrowLength = 0.2,
                alpha = 0.2,
                axis1 = axis1,
                axis2 = axis2) +
      addCircleOfCor() + xlab(paste("Component", axis1)) + ylab(paste("Component", axis2))

   print(LoadingMapWithCircles)

   signed.ctrJ <- resPCA$Fixed.Data$ExPosition.Data$cj *
      sign(resPCA$Fixed.Data$ExPosition.Data$fj)

   # plot contributions for component 1
   ctrJ.1 <- PrettyBarPlot2(signed.ctrJ[,axis1],
                            threshold = 1 / NROW(signed.ctrJ),
                            font.size = 3,
                            #color4bar = gplots::col2hex(res_pcaInf$Fixed.Data$Plotting.Data$fj.col),
                            ylab = 'Contributions',
                            ylim = c(1.2*min(signed.ctrJ[,axis1:axis2]), 1.2*max(signed.ctrJ[,axis1:axis2]))) +
      ggtitle("Contribution barplots", subtitle = paste("Component", axis1))

   # plot contributions for component 2
   ctrJ.2 <- PrettyBarPlot2(signed.ctrJ[,axis2],
                            threshold = 1 / NROW(signed.ctrJ),
                            font.size = 3,
                            #color4bar = gplots::col2hex(res_pcaInf$Fixed.Data$Plotting.Data$fj.col),
                            ylab = 'Contributions',
                            ylim = c(1.2*min(signed.ctrJ[,axis1:axis2]), 1.2*max(signed.ctrJ[,axis1:axis2]))) +
      ggtitle("",subtitle = paste("Component", axis2))

   print(ctrJ.1)
   print(ctrJ.2)

}

#' PCA with inference and graphs.
#' @param data The numerical data for PCA.
#' @param design (default = NULL) A vector or dummy-coded matrix that describes the rows.
#' @param make_design_nominal (default = TRUE) If design is a dummy-coded matrix, should be set to FALSE.
#' @param col4obs (default = "olivedrab3") A single color or vector of colors whose length is equal to nrow(data).
#' @param col4group (default = "olivedrab3") A single color or vector of colors whose length is the number of groups in design.
#' @param center (default = TRUE) Whether to center variables
#' @param scale (default = "SS1") Whether to scale variables
#' @param want34 (default = FALSE) By default, only prints dimensions 1 and 2. Set to TRUE for 3 and 4.
#' @param inference (default = TRUE) When design contains 2 or more groups, computes CI and TI on observations. FALSE if no groups.
#' @importFrom corrplot corrplot
#' @importFrom PTCA4CATA PlotScree
#' @import InPosition
#' @export
mora_pca <- function(data,
                    design = NULL,
                    make_design_nominal = TRUE,
                    col4obs = "olivedrab3",
                    col4group = "olivedrab3",
                    center = TRUE,
                    scale = "SS1",
                    want34 = FALSE,
                    inference = TRUE){
   data_cor <- cor(data)

   corrplot(data_cor, tl.cex = 0.7, tl.pos = "lt", tl.col = "black",
            addCoefasPercent = TRUE, addCoef.col = "black",
            number.cex = 0.5, method = "color")

   resPCA <- epPCA.inference.battery(data, center = center, scale = scale,
                                     DESIGN = design, graphs = FALSE, test.iters = 1000)

   scree <- PlotScree(ev = resPCA$Fixed.Data$ExPosition.Data$eigs,
                           p.ev = resPCA$Inference.Data$components$p.vals,
                           plotKaiser = TRUE)

   pca_fscores(resPCA = resPCA, design = design, axis1 = 1, axis2 = 2,
               col4obs = col4obs, col4group = col4group, inference = inference)

   pca_columns(resPCA, data = data, 1, 2)

   if(want34){
      pca_fscores(resPCA = resPCA, design = design, axis1 = 3, axis2 = 4,
                  col4obs = col4obs, col4group = col4group, inference = inference)

      pca_columns(resPCA, data = data, 3, 4)
   }
}
