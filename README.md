# Data Analysis using the Morse-Smale Complex
- Paper reference: Chen, Yen-Chi, Christopher R. Genovese, and Larry Wasserman. "Statistical Inference using the Morse-Smale Complex." arXiv preprint arXiv:1506.08826 (2015).
- `MSHD.R`: the main script for the method.
- `example.R`: an example (used in the paper) for applying MSHD.

## MSHD.R

### MSHD
`MSHD = function(grids, grids.values, pLevel=0.01, knn=3*ncol(grids)) UseMethod("MSHD")`

`MSHD.default = function(grids, grids.values, pLevel=0.01, knn=3*ncol(grids))`
- Visualization for a high dimensional function
- Inputs:
  - grids: The grid points where the high dimensional function is evaluated.
  - grids.values: The values for the high dimensional function at the corresponding grid point.
  - pLevel: The persistent level used in $MSR$ package.
  - knn: The number of nearest neighbor used in $MSR$ package.
- Output:
  - A S4 object of class "MSHD". A list consisting:
    - grids: The input grid points.
    - grids.values: The input grid values.
    - label.cell: The cell label for each grid point.
    - cell.info: A matrix containing information for each cell. Each row corresponds to a cell.
    - beta: The regression coefficients for linear approximation for the high dimensional function.

### cell.plot
`cell.plot = function(x, lv.pos = NULL, lv.neg = NULL, grid.width = NULL, r_circle = 1, col.cell = "limegreen", col.cell.filled = "limegreen", col.edge = "brown", col.pie = c("dodgerblue", "red"), lwd.cell = 4, lwd.edge = 10, use.legend=T)`

#' Visualizing a high dimensional function using Morse-Smale cells.
#' @param x An MSHD object.
#' @param lv.pos The threshold for claiming a grid point is significantly high.
#' @param lv.neg The threshold for claiming a grid point is significantly low.
#' @param grid.width The grid width that will be used to compute the number of grid
#' points that two cells share.
#' @param r_circle The visualization radius for each cell.
#' @param col.cell Color for the cell (boundaries of the circle).
#' @param col.cell.filled The filled-in color for the cell.
#' @param col.edge Color for the edge.
#' @param col.pie Color for significantly positive and negative region.
#' @param lwd.cell Width of the circle for the cell.
#' @param lwd.edge Width for the edges.
#' @param cex.cell.txt `txt' for the text label for each cell.
#' @param pos.cell.txt `pos' for the text label for each cell.
#' @param offset.cell.txt `offset' for the text label for each cell.
#' @param use.legend To display the legend or not.
#' 
#' @return A plot for visualization using cells and a list consisting 
#' \item{cell.vis}
#' The visualization coordinates for each cell.
#' \item{cell.edgepoints}
#' Number of edge points for a pair of cell.
#' \item{cell.size}
#' The size (number of grid points) for each cell, including total size,
#' significantly positive grid points, and significantly negative grid points.
#' 
#' @export
