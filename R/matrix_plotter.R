#' Different default methods to visualize matrices
#' @param data The matrix to be visualized.
#' @param method Current options: 'corrplot', 'heatmap'.
#' @param col (default = NULL) Colors for the values to be plotted.
#' @param coef (default = TRUE) Whether to display the values of the matrix
#' @param is.corr (default = TRUE) Whether the matrix is a correlation matrix
#' @export
matrix_plotter <- function(data, method, col = NULL, coef = TRUE, is.corr = TRUE ){

   method.types <- c("corrplot", "heatmap")

   if(!(method %in% method.types)){
      warning(paste0("method ", method, " not found. Setting method to corrplot."))
      method <- "corrplot"
   }

   n.row <- NROW(data)
   n.col <- NCOL(data)

   mindim <- min(n.row, n.col)
   maxdim <- max(n.row, n.col)

   if(coef & maxdim > 60){
      message("Printing coefficients not recommended for nvar > 60, setting coef = FALSE")
      coef <- FALSE
   }

   if(is.null(col)){
      col2 <- grDevices::colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                 "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                 "#4393C3", "#2166AC", "#053061"))
      col <- col2(200)
   }


   if(method == "corrplot"){
      corrplot::corrplot(data,
                         tl.cex = 4 * (1/sqrt(maxdim)),
                         tl.pos = "lt",
                         tl.col = "black",
                         addCoefasPercent = TRUE,
                         addCoef.col = if(coef) "black" else NULL,
                         number.cex = 4 * (1/sqrt(maxdim)),
                         method = "color",
                         is.corr = is.corr,
                         col = col)
   }
   else if(method == "heatmap"){
      if(coef){
         gplots::heatmap.2(data,
                           Rowv = FALSE,
                           Colv = FALSE,
                           dendrogram = "none",
                           col = col,
                           trace = "none",
                           key = TRUE,
                           margins = c(1,1),
                           lmat = rbind(c(4,3), c(1,2)),
                           lwid = c(5, 0.1),
                           lhei = c(0.5, 5),
                           cexRow = 1/log10(n.row),
                           cexCol = 1/log10(n.col),
                           cellnote = if(is.corr) round(data*100) else signif(data, 2),
                           notecol = "black",
                           notecex = 4 * (1/sqrt(maxdim)),
                           key.title = NA,
                           key.xlab = NA,
                           density.info = 'none'
                           )
      }
      else{
         gplots::heatmap.2(data,
                           Rowv = FALSE,
                           Colv = FALSE,
                           dendrogram = "none",
                           col = col,
                           trace = "none",
                           key = TRUE,
                           margins = c(1,1),
                           lmat = rbind(c(4,3), c(1,2)),
                           lwid = c(5, 0.1),
                           lhei = c(0.5, 5),
                           cexRow = 1/log10(n.row),
                           cexCol = 1/log10(n.col),
                           key.title = NA,
                           key.xlab = NA,
                           density.info = 'none'
                           )
      }
   }
}
