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
- Visualizing a high dimensional function using Morse-Smale cells.
- Inputs:
  - x: An MSHD object.
  - lv.pos: The threshold for claiming a grid point is significantly high.
  - lv.neg: The threshold for claiming a grid point is significantly low.
  - grid.width: The grid width that will be used to compute the number of grid points that two cells share.
  - r_circle: The visualization radius for each cell.
  - col.cell: Color for the cell (boundaries of the circle).
  - col.cell.filled: The filled-in color for the cell.
  - col.edge: Color for the edge.
  - col.pie: Color for significantly positive and negative region.
  - lwd.cell: Width of the circle for the cell.
  - lwd.edge: Width for the edges.
  - cex.cell.txt: `txt` of the text label for each cell.
  - pos.cell.txt: `pos` of the text label for each cell.
  - offset.cell.txt: `offset` of the text label for each cell.
  - use.legend: To display the legend or not.
- Output:
  - A plot for visualization using cells and a list consisting:
    - cell.vis: The visualization coordinates for each cell.
    - cell.edgepoints: Number of edge points for a pair of cell.
    - cell.size: The size (number of grid points) for each cell, including total size, significantly positive grid points, and significantly negative grid points.


### MSE.test
`MSEtest = function(data1, data2, grids, h=NULL, pLevel=0.01, knn=3*ncol(grids))`
- Morse-Smale Two Sample Test using Energy Statistics and KDE.
- Inputs:
  - data1: The first group of data
  - data2: The second group of data
  - grids: The grid points where the density will be evaluated
  - h: The smoothing parameter.
  - pLevel: The persistent level used in `MSR` package.
  - knn: The number of nearest neighbor used in `MSR` package.
- Output:
  - A S4 object of class "MSHD". A list consisting:
    - cell1.data: The sample from 1st group for constructing Morse-Smale cells.
    - cell2.data: The sample from 2nd group for constructing Morse-Smale cells.
    - test1.data: The sample from 1st group for conducting energy test.
    - test2.data: The sample from 2nd group for conducting energy test.
    - grids: The grid points.
    - grids.density: A matrix containing the density evaluated at the grid points. First column is the density from cell1 data; second column is from cell2 data; and the last column is the density difference.
    - labels: A list consisting of cell labels for grids, cell1, cell2, test1, and test2.
    - p.values: The p.values for each cell. The value is \emph{NA} if there is only one point from either test1 or test2 for this cell.
    - bandwidth: The smoothing bandwidth used to estimate density.


