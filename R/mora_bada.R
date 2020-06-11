#' BADA with inference and graphs.
#' @param data Numerical data table with observations on rows and variables on columns.
#' @param design A vector or dummy-coded matrix indicating group membership of rows.
#' @param col4obs A single color or vector of colors whose length is equal to nrow(data).
#' @param col4group A single color or vector of colors whose length is the number of groups in design.
#' @param center (default = TRUE) Whether to center variables
#' @param scale (default = "SS1") Whether to scale variables
#' @param make_design_nominal (default = TRUE) If design is a dummy-coded matrix, should be set to FALSE.
#' @param test.iters (default = 1000) Number of iterations for bootstrap.
#' @param important (default = FALSE) If TRUE, graphs have only the important saliences/bootstrap ratios.
#' @export
mora_bada <- function(data,
                      design,
                      col4obs,
                      col4group,
                      center = TRUE,
                      scale = "SS1",
                      make_design_nominal = TRUE,
                      test.iters = 1000,
                      important = FALSE){

   data.scaled <- ExPosition::expo.scale(data, center, scale)
   data.means <- PTCA4CATA::getMeans(data.scaled, design)

   # TO DO: make this look better (different color scheme)
   stats::heatmap(as.matrix(data.means), Rowv = NA, Colv = NA, scale = "column", cexCol = 0.5)


   resBADA.inf <- TInPosition::tepBADA.inference.battery(data, center = center,
                                            scale = scale, DESIGN = design,
                                            graphs = FALSE,
                                            make_design_nominal = make_design_nominal,
                                            test.iters = test.iters)

   myScree <- PTCA4CATA::PlotScree(ev = resBADA.inf$Fixed.Data$TExPosition.Data$eigs,
                        p.ev = resBADA.inf$Inference.Data$components$p.vals,
                        plotKaiser = TRUE)

   Obsmap <- PTCA4CATA::createFactorMap(
      resBADA.inf$Fixed.Data$TExPosition.Data$fii,
      col.points = col4obs,
      col.labels = col4obs,
      display.labels = FALSE,
      alpha.points = .5
   )


   #rownames(resBADA.inf$Fixed.Data$TExPosition.Data$fi) <- removeFirstChar(rownames(resBADA.inf$Fixed.Data$TExPosition.Data$fi))

   fi.labels <- PTCA4CATA::createxyLabels.gen(lambda = resBADA.inf$Fixed.Data$TExPosition.Data$eigs,
                                   tau = round(resBADA.inf$Fixed.Data$TExPosition.Data$t),
                                   axisName = "Component ")

   fii.map <- Obsmap$zeMap + fi.labels
   #print(fii.map)

   #fii.means <- getMeans(resBADA.inf$Fixed.Data$TExPosition.Data$fii, design)
   Fi <- resBADA.inf$Fixed.Data$TExPosition.Data$fi

   groupMap <- PTCA4CATA::createFactorMap(Fi,
                               # use the constraint from the main map
                               constraints = Obsmap$constraints,
                               col.points = col4group[rownames(Fi)],
                               cex = 5,  # size of the dot (bigger)
                               col.labels = col4group[rownames(Fi)],
                               text.cex = 6,
                               alpha.points = 0.8)

   groupAndObsMap <- Obsmap$zeMap + fi.labels + groupMap$zeMap_dots + groupMap$zeMap_text

   #print(groupAndObsMap)

   #rownames(resBADA.inf$Inference.Data$boot.data$fi.boot.data$boots) <- removeFirstChar(rownames(resBADA.inf$Inference.Data$boot.data$fi.boot.data$boots[,c(1,2),]))

   CIEllip <- PTCA4CATA::MakeCIEllipses(resBADA.inf$Inference.Data$boot.data$fi.boot.data$boots[,c(1,2),],
                             col = col4group[rownames(resBADA.inf$Inference.Data$boot.data$fi.boot.data$boots[,c(1,2),])])


   FiWithCI <- groupAndObsMap + CIEllip

   print(FiWithCI)

   TIplot <- PTCA4CATA::MakeToleranceIntervals(resBADA.inf$Fixed.Data$TExPosition.Data$fii,
                                    design = design,
                                    names.of.factors =  c("Dim1","Dim2"), # needed
                                    col = col4group,
                                    line.size = .50,
                                    line.type = 3,
                                    alpha.ellipse = .2,
                                    alpha.line    = .4,
                                    p.level       = .95)

   fiWithTI <- groupAndObsMap + TIplot

   print(fiWithTI)

   BR <- resBADA.inf$Inference.Data$boot.data$fj.boot.data$tests$boot.ratios

   impBR1 <- vector()
   if(important){
      impVar1 <- which(resBADA.inf$Fixed.Data$TExPosition.Data$cj[,1] >
                          (1/nrow(resBADA.inf$Fixed.Data$TExPosition.Data$cj)))
      impVar2 <- which(resBADA.inf$Fixed.Data$TExPosition.Data$cj[,2] >
                          (1/nrow(resBADA.inf$Fixed.Data$TExPosition.Data$cj)))
      impVars <- union(impVar1, impVar2)

      fjMap <- PTCA4CATA::createFactorMap(resBADA.inf$Fixed.Data$TExPosition.Data$fj[impVars,],
                               pch = 19,
                               cex = 2,
                               text.cex = 2.5,
                               alpha.points = 0.7,
                               display.labels = TRUE
                               #col.points = resBADA.inf$Fixed.Data$Plotting.Data$fj.col,
                               #col.labels = resBADA.inf$Fixed.Data$Plotting.Data$fj.col
      )

      arrows <- PTCA4CATA::addArrows(resBADA.inf$Fixed.Data$TExPosition.Data$fj[impVars,],
                          color = "black", alpha = 0.3)


      impBR1 <- which(abs(BR[,1]) > 2)
      impBR2 <- which(abs(BR[,2]) > 2)

      ba001.BR1 <- PTCA4CATA::PrettyBarPlot2(BR[impBR1,1],
                                  threshold = 2,
                                  font.size = 3,
                                  #color4bar = gplots::col2hex(resBADA.inf$Fixed.Data$Plotting.Data$fj.col),
                                  ylab = 'Bootstrap ratios') +
         ggtitle("Bootstrap ratios", subtitle = "Component 1")

      # Plot the bootstrap ratios for Dimension 2

      ba002.BR2 <- PTCA4CATA::PrettyBarPlot2(BR[impBR2,2],
                                  threshold = 2,
                                  font.size = 3,
                                  #color4bar = gplots::col2hex(resBADA.inf$Fixed.Data$Plotting.Data$fj.col),
                                  ylab = 'Bootstrap ratios') +
         ggtitle("Bootstrap ratios",subtitle = "Component 2")

   }
   else{
      fjMap <- PTCA4CATA::createFactorMap(resBADA.inf$Fixed.Data$TExPosition.Data$fj,
                               pch = 19,
                               cex = 2,
                               text.cex = 2.5,
                               alpha.points = 0.7,
                               display.labels = TRUE
                               #col.points = resBADA.inf$Fixed.Data$Plotting.Data$fj.col,
                               #col.labels = resBADA.inf$Fixed.Data$Plotting.Data$fj.col
      )

      arrows <- PTCA4CATA::addArrows(resBADA.inf$Fixed.Data$TExPosition.Data$fj, color = "black", alpha = 0.3)

      ba001.BR1 <- PTCA4CATA::PrettyBarPlot2(BR[,1],
                                  threshold = 2,
                                  font.size = 3,
                                  #color4bar = gplots::col2hex(resBADA.inf$Fixed.Data$Plotting.Data$fj.col),
                                  ylab = 'Bootstrap ratios') +
         ggtitle("Bootstrap ratios", subtitle = "Component 1")

      # Plot the bootstrap ratios for Dimension 2

      ba002.BR2 <- PTCA4CATA::PrettyBarPlot2(BR[,2],
                                  threshold = 2,
                                  font.size = 3,
                                  #color4bar = gplots::col2hex(resBADA.inf$Fixed.Data$Plotting.Data$fj.col),
                                  ylab = 'Bootstrap ratios') +
         ggtitle("Bootstrap ratios", subtitle = "Component 2")

   }

   fj.plot <- fjMap$zeMap + fi.labels + arrows

   print(fj.plot)

   print(ba001.BR1)
   print(ba002.BR2)

   print(glue::glue(" \n \nOmnibus p val = {resBADA.inf$Inference.Data$omni$p.val}.\n",
        "R squared = {round(resBADA.inf$Inference.Data$r2$r2, 2)}, p = {resBADA.inf$Inference.Data$r2$p.val}.\n",
        "Fixed effect accuracy = {round(resBADA.inf$Inference.Data$loo.data$fixed.acc, 2)}. Confusion matrix:\n"))
   print(resBADA.inf$Inference.Data$loo.data$fixed.confuse)
   print(glue::glue("Random effect accuracy = {round(resBADA.inf$Inference.Data$loo.data$loo.acc, 2)}. Confusion matrix:\n"))
   print(resBADA.inf$Inference.Data$loo.data$loo.confuse)

   # TO DO: fix for a consistent return value
   if(important){
      return(impBR1)
   }
   #return(resBADA.inf)
}
