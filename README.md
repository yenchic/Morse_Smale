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

###

