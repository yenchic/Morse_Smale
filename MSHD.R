# Author: Yen-Chi Chen, Christopher R. Genovese, Larry Wasserman
# Maintainer: Yen-Chi Chen <ga014528@gmail.com>
# Reference: Statistical Inference using Morse-Smale Complex
# Date: 06/22/2015
#' @import msr
library(msr)
#' @import mapplots
library(mapplots)
#' @import plotrix
library(plotrix)
#' @import RANN
library(RANN)
#' @import TDA
library(TDA)
#' @import energy
library(energy)

#' Visualization for a high dimensional function.
#' @param grids The grid points where the high dimensional function is evaluated.
#' @param grids.values The values for the high dimensional function 
#' at the corresponding grid point.
#' @param pLevel The persistent level used in `MSR' package.
#' @param knn The number of nearest neighbor used in `MSR' package.
#' @export
MSHD = function(grids, grids.values, pLevel=0.01, knn=3*ncol(grids)) UseMethod("MSHD")

#' Visualization for a high dimensional function
#' @param grids The grid points where the high dimensional function is evaluated.
#' @param grids.values The values for the high dimensional function 
#' at the corresponding grid point.
#' @param pLevel The persistent level used in `MSR' package.
#' @param knn The number of nearest neighbor used in `MSR' package.
#' 
#' @return A S4 object of class "MSHD". A list consisting
#' \item{grids} 
#' The input grid points.
#' \item{grids.values}
#' The input grid values.
#' \item{label.cell}
#' The cell label for each grid point.
#' \item{cell.info}
#' A matrix containing information for each cell. Each row corresponds to a cell.
#' \item{beta}
#' The regression coefficients for linear approximation for the 
#' high dimensional function.
#' @export
MSHD.default = function(grids, grids.values, pLevel=0.01, knn=3*ncol(grids)){
  d = ncol(grids)
  result = list()
  
  if(nrow(grids)!=length(grids.values)){
    cat("Number of grid points does not match number of grid values!")
    stop
  }
  ms0 = msc.nn(x=grids, y=grids.values, pLevel=pLevel, knn=knn)
  
  idx_cl = ms0$level[[1]]$partition
    ### Labels for cell
  
  modes_idx = unique(ms0$level[[1]]$max)
  mini_idx = unique(ms0$level[[1]]$min)
  modes = grids[modes_idx,]
  minis = grids[mini_idx,]
    ### Local modes and minima
  
    ### The piecewise linear regression
  cell_beta= matrix(NA, nrow=max(idx_cl), ncol=d+1)
  for(i in 1:max(idx_cl)){
    idx_tmp = which(idx_cl==i)
    
    fit_tmp = lm(grids.values[idx_tmp]~as.matrix(grids[idx_tmp,]))
    cell_beta[i,] = fit_tmp$coeff
  }
  
  result$grids = grids
  result$grids.values = grids.values
  
  result$label.cell = idx_cl
  
  result$cell.info = cbind(ms0$level[[1]]$max, 
                           ms0$level[[1]]$min, 
                           ms0$level[[1]]$partitionSize)
  rownames(result$cell.info) = paste("cell_",1:max(idx_cl), sep="")
  colnames(result$cell.info) = c("Index_modes","Index_minima", "Size")
  
  result$beta = cell_beta
  rownames(result$beta) = paste("cell_", 1:max(idx_cl), sep="")
  colnames(result$beta) = c("Intercept", colnames(grids))

  class(result) = "MSHD"
  return(result)
}
print.MSHD = function(x){
  cat("Cells information:\n")
  print(x$cell.info)
  cat("Linear approximation coefficients:\n")
  print(x$beta)
}
                     
#' Plotting a high dimensional function using Morse-Smale signature.
#' @param x An MSHD object.
#' @param col.modes The color for local modes (nodes).
#' @param col.minima The color for local minima (nodes).
#' @param col.cells The color for cells (edges).
#' @param cex.modes `cex' for local modes.
#' @param cex.minima `cex' for local minima.
#' @param pch.modes `pch' for local modes.
#' @param pch.minima `pch' for local minima.
#' @param lwd.edge `lwd' for edges. This is for the edge (cell) with highest 
#' slope (L2 norm of linear approximation coefficients). The other edges will
#' have thiner width.
#' @param use.legend To display the legend or not.
#' 
#' @return A plot for visualization and a list consisting
#' \item{modes.vis}
#' The low dimensional representation for the position of each local mode.
#' \item{minima.vis}
#' The low dimensional representation for the position of each local minimum.
#' \item{edges}
#' To show what are the nodes (pair of local mode and minimum) that each edge 
#' (cell) connect.
#' \item{beta_slope}
#' The L2 norm (without intercept) for the regression coefficients for each cell.
#' 
#' @export
plot.MSHD = function(x, col.modes="dodgerblue", col.minima = "limegreen", 
                     col.cells = "brown", cex.modes = 12, cex.minima = 8,
                     pch.modes = 20, pch.minima = 20, lwd.edge = 25, use.legend=T, 
                     ...){
  
    ### points (critical points)
  modes_idx = unique(x$cell.info[,1])
  minima_idx = unique(x$cell.info[,2])
  modes = x$grids[modes_idx,]
  minima = x$grids[minima_idx,]
  CP = rbind(modes, minima)
  CP_mds = cmdscale(dist(CP))
  modes_mds = CP_mds[1:nrow(modes),]
  minima_mds = CP_mds[(nrow(modes)+1):(nrow(CP)),]
  
  
    ### edges (cells)
  cell_beta_s = rowSums(x$beta[,-1]^2)
  cell_beta_s = sqrt(cell_beta_s/max(cell_beta_s))
  edges = matrix(NA, nrow=2, ncol=nrow(x$cell.info))
  for(i in 1:nrow(x$cell.info)){
    for(j in 1:nrow(x$cell.info)){
      a = which(modes_idx==x$cell.info[i,1])
      b = which(minima_idx==x$cell.info[j,2])
      edges[1,i] = a
      edges[2,j] = b
    }
  }
  
  x.width = max(CP_mds[,1]) - min(CP_mds[,1])
  y.width = max(CP_mds[,2]) - min(CP_mds[,2])
  
  par(mar=c(0.5,0.5,0.5,0.5))
  plot(CP_mds, xlim = c(min(CP_mds[,1])-0.2*x.width,max(CP_mds[,1])+0.2*x.width),
       ylim = c(min(CP_mds[,2])-0.2*y.width, max(CP_mds[,2])+0.2*y.width), 
       xlab="", ylab="", xaxt="n",yaxt="n")
  for(i in 1:ncol(edges)){
    segments(x0=modes_mds[edges[1,i],1], y0=modes_mds[edges[1,i],2], 
             x1=minima_mds[edges[2,i],1], y1=minima_mds[edges[2,i],2], 
             lwd = lwd.edge*cell_beta_s[i]+0.5, col="brown")
  }
  points(modes_mds, col=col.modes, pch=pch.modes, cex=cex.modes)
  points(minima_mds, col=col.minima, pch=pch.minima, cex=cex.minima)
  if(use.legend){
    legend("bottomleft",c("Modes","Minima"), col=c("dodgerblue","limegreen"), 
         pch=c(pch.modes, pch.minima), pt.cex=c(cex.modes,cex.minima)*0.4, cex=1)
    legend("topleft", legend="Cells", col=col.cells, lwd=lwd.edge*0.6, cex=1)
  }
  
  result = list()
  result$modes.vis = modes_mds
  rownames(result$modes.vis) = paste("Mode_", 1:length(modes_idx), sep="")
  
  result$minima.vis = minima_mds
  rownames(result$minima.vis) = paste("Minimum_", 1:length(minima_idx), sep="")
  
  result$edges = edges
  rownames(result$edges) = c("Mode_label", "Minimum_label")
  colnames(result$edges) = paste("cell_", 1:nrow(x$cell.info), sep="")
  
  result$beta_slope = cell_beta_s
  
  class(result)= "MSHD.plot"
  return(result)
}
print.MSHD.plot = function(x){
  cat("")
}


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
cell.plot = function(x, lv.pos = NULL, lv.neg = NULL, grid.width = NULL, 
                     r_circle = 1, 
                     col.cell = "limegreen", col.cell.filled = "limegreen",
                     col.edge = "brown", 
                     col.pie = c("dodgerblue", "red"), 
                     lwd.cell = 4, lwd.edge = 10, 
                     use.legend=T
                     ){
  if(class(x)!="MSHD"){
    cat("x must be a MSHD object!")
    stop
  }
    ### MDS for the mean position of each cell
  cell_c = NULL
  for(idx in 1:nrow(x$cell.info)){
    idx_tmp = which(x$label.cell==idx)
    
    cell_tmp = colSums(x$grids[idx_tmp,])/length(idx_tmp)
    cell_c = rbind(cell_c,cell_tmp)
  }
  cell_mds = cmdscale(dist(cell_c))
  cell_size_s = sqrt(x$cell.info[,3]/max(x$cell.info[,3]))
  
    ### significance regions
  cell_size_psig = rep(0, nrow(x$cell.info))
  cell_size_nsig = rep(0, nrow(x$cell.info))
  if(!is.null(lv.pos)){
    idx_psig = which(x$grids.values > lv.pos)
    idx_cl_psig = x$label.cell[idx_psig]
    for(i in 1:nrow(x$cell.info)){
      cell_size_psig[i] = length(which(idx_cl_psig==i))
    }
  }
  if(!is.null(lv.neg)){
    idx_nsig = which(x$grids.values < lv.neg)
    idx_cl_nsig = x$label.cell[idx_nsig]
    for(i in 1:nrow(x$cell.info)){
      cell_size_nsig[i] = length(which(idx_cl_nsig==i))
    }
  }
  cell_r_psig = cell_size_psig/x$cell.info[,3]
  cell_r_nsig = cell_size_nsig/x$cell.info[,3]
  
  
    ### edges
  if(is.null(grid.width)){
    grid.width = sqrt(ncol(x$grids))*sqrt(sum((x$grids[1,]-x$grids[2,])^2))
  }
  cell_list = list()
  cell_edge = matrix(0, nrow=nrow(x$cell.info), ncol=nrow(x$cell.info))
  for(idx in 1:nrow(x$cell.info)){
    idx_tmp = which(x$label.cell==idx)
    cell_list[[idx]]= x$grids[idx_tmp,]
  }
  
  for(idx1 in 1:nrow(x$cell.info)){
    for(idx2 in 1:nrow(x$cell.info)){
      cell_A = cell_list[[idx1]]
      cell_B = cell_list[[idx2]]
      nn_dist = nn2(cell_A, cell_B, k=1)$nn.dist
      cell_edge[idx1,idx2] = length(which(nn_dist<grid.width))
    }
  }
  cell_edge_s = cell_edge/x$cell.info[,3]
  cell_edge_ss = (cell_edge_s+t(cell_edge_s))/2

  
    ### making plot
  x.width = max(cell_mds[,1]) - min(cell_mds[,1])
  y.width = max(cell_mds[,2]) - min(cell_mds[,2])
  
  par(mar=c(0.5,0.5,0.5,0.5))
  plot(cell_mds, xlim = c(min(cell_mds[,1])-0.2*x.width, 
                          max(cell_mds[,1])+0.2*x.width),
       ylim = c(min(cell_mds[,2])-0.2*y.width, max(cell_mds[,2])+0.2*y.width), 
       xlab="", ylab="", xaxt="n",yaxt="n", col="white")
    # modes (edges)
  for(idx1 in 1:(nrow(x$cell.info)-1)){
    for(idx2 in (idx1+1):nrow(x$cell.info)){
      if(cell_edge_ss[idx1,idx2]>0.1){
        segments(x0= cell_mds[idx1,1], y0= cell_mds[idx1,2], x1= cell_mds[idx2,1], 
                 y1= cell_mds[idx2, 2], lwd= lwd.edge*cell_edge_ss[idx1,idx2], 
                 col=col.edge)
      }
    }
  }
    # regions (cells)
  min_cell = min(dist(cell_mds))
  for(i in 1:nrow(x$cell.info)){
    draw.circle(x=cell_mds[i,1], y=cell_mds[i,2], 
                radius=r_circle*min_cell*cell_size_s[i]/2,
                lwd = lwd.cell, border=col.cell, col=col.cell.filled)
    if(!(is.null(lv.pos)&is.null(lv.neg))){
      r_sig = cell_r_psig[i]+ cell_r_nsig[i]
      if(r_sig>0){
        add.pie(c(cell_r_psig[i]/r_sig, cell_r_nsig[i]/r_sig), 
                radius=r_circle*min_cell*r_sig*cell_size_s[i]/2, 
                x=cell_mds[i,1], y=cell_mds[i,2], col=col.pie, 
                labels=c("",""))
      }
    }
  }
  if(use.legend&!(is.null(lv.pos)&is.null(lv.neg))){
    legend("bottomright",legend=c("Sig. Positive", "Sig. Negative"), 
           col=col.pie, cex=1, pch=20, pt.cex=4)
  }
  
  ### output
  result = list()
  result$cell.vis = cell_mds
  rownames(result$cell.vis) = paste("cell_", 1:nrow(x$cell.info), sep="")
  
  result$cell.edgepoints = cell_edge
  rownames(result$cell.edgepoints) = paste("cell_", 1:nrow(x$cell.info), sep="")
  colnames(result$cell.edgepoints) = paste("cell_", 1:nrow(x$cell.info), sep="")
  
  result$cell.size = cbind(x$cell.info[,3], cell_r_psig, cell_r_nsig)
  colnames(result$cell.size) = c("Size", "Positive","Negative")
  
  class(result) = "MSHD.cell.plot"
  return(result)
}
print.MSHD.cell.plot = function(x){
  cat("")
}


#' Morse-Smale Two Sample Test using Energy Statistics and KDE.
#' @param data1 The first group of data
#' @param data2 The second group of data
#' @param grids The grid points where the density will be evaluated
#' @param h The smoothing parameter.
#' @param pLevel The persistent level used in `MSR' package.
#' @param knn The number of nearest neighbor used in `MSR' package.
#' 
#' @return A S4 object of class "MSHD". A list consisting
#' \item{cell1.data} 
#' The sample from 1st group for constructing Morse-Smale cells.
#' \item{cell2.data} 
#' The sample from 2nd group for constructing Morse-Smale cells.
#' \item{test1.data} 
#' The sample from 1st group for conducting energy test.
#' \item{test2.data} 
#' The sample from 2nd group for conducting energy test.
#' \item{grids}
#' The grid points.
#' \item{grids.density}
#' A matrix containing the density evaluated at the grid points.
#' First column is the density from cell1 data; second column
#' is from cell2 data; and the last column is the density difference.
#' \item{labels}
#' A list consisting of cell labels for grids, cell1, cell2, test1, and test2.
#' \item{p.values}
#' The p.values for each cell. The value is \emph{NA} if there is only one 
#' point from either test1 or test2 for this cell.
#' \item{bandwidth}
#' The smoothing bandwidth used to estimate density.
#' @export
MSEtest = function(data1, data2, grids, h=NULL, pLevel=0.01, knn=3*ncol(grids)){
  n_1 = nrow(data1)
  n_2 = nrow(data2)
  n_all = (n_1+n_2)/2
  d = ncol(data1)
  
  idx_11 = sample(n_1, floor(n_1/2))
  idx_12 = c(1:n_1)[-idx_11]
  idx_21 = sample(n_2, floor(n_2/2))
  idx_22 = c(1:n_2)[-idx_21]
  
    # for cell
  C1 = data1[idx_11,]
  C2 = data2[idx_21,]
  C_all = rbind(C1,C2)
    
    # for test
  T1 = data1[idx_12,]
  T2 = data2[idx_22,]

  ### KDE for sample 1
  if(is.null(h)){
    h = (4/(d+4))^(1/(d+6))/n_all^(1/(d+6))*mean(apply(C_all,2,sd))
  }
  grids.C1den = kde(C1, grids, h)
  grids.C2den = kde(C2, grids, h)
  
  grids.diff = grids.C1den-grids.C2den
    ### Building cells (for grids)
  ms0 = msc.nn(x=grids, y=grids.diff/max(abs(grids.diff)), pLevel=pLevel, knn=knn)
  
    ### Assign cell labels to C1 and C2
  cl_C1 = ms0$level[[1]]$partition[nn2(data= grids, query = C1,k=1)$nn.idx]
  cl_C2 = ms0$level[[1]]$partition[nn2(data= grids, query = C2,k=1)$nn.idx]
  
    ### Assign cell labels to T1 and T2
  cl_T1 = ms0$level[[1]]$partition[nn2(data= grids, query = T1,k=1)$nn.idx]
  cl_T2 = ms0$level[[1]]$partition[nn2(data= grids, query = T2,k=1)$nn.idx]
  
  
    ### Energy Test
  p_m = rep(NA, max(ms0$level[[1]]$partition))
  for(i_e in 1:max(ms0$level[[1]]$partition)){
    T1_tmp = T1[which(cl_T1==i_e),]
    T2_tmp = T2[which(cl_T2==i_e),]	
    if(length(T1_tmp)>d&length(T2_tmp)>d){
      T_pull = rbind(T1_tmp,T2_tmp)
      E_test = eqdist.etest(T_pull, c(nrow(T1_tmp), nrow(T2_tmp)))
      p_m[i_e] = E_test$p.value
    }
  }
  result = list()
  
  result$cell1.data = C1
  result$cell2.data = C2
  result$test1.data = T1
  result$test2.data = T2
  result$grids = grids
  result$grids.density = cbind(grids.C1den, grids.C2den, grids.diff)
  colnames(result$grids.density) = c("C1.density","C2.density","diff.")
  result$labels = list()
  result$labels$grids = ms0$level[[1]]$partition
  result$labels$cell1 = cl_C1
  result$labels$cell2 = cl_C2
  result$labels$test1 = cl_T1
  result$labels$test2 = cl_T2
  result$p.values = p_m
  result$bandwidth = h
  
  class(result) = "MSE.test"
  return(result)
}
print.MSE.test = function(x){
  cat("Size for sample 1:\n")
  print(nrow(x$cell1.data)+nrow(x$test1.data))
  cat("Size for sample 2:\n")
  print(nrow(x$cell2.data)+nrow(x$test2.data))
  cat("P-values for each cell:\n")
  print(x$p.values)
}
  
#' Four Gaussian mixture
#' @param n Sample size to generate.
#' @param par.pi The proportion of each mixture. A vector consists of 4 elements.
#' Note that the centers for the 4 mixture is (0,0), (1,0), (0,1) and (1,1).
#' @param par.err The standard deviation of each mixture. A vector consists of 4 elements.
#' @return A n*2 matrix consists of the 4 Gaussian mixture.
#' @export
GMM4 = function(n=1000, par.pi=c(0.2,0.5,0.2,0.1), par.err = c(0.2,0.2,0.2,0.2)){
  X1 = cbind(rnorm(n*par.pi[1], 0, sd=par.err[1]), rnorm(n*par.pi[1], 0, sd=par.err[1]))
  X2 = cbind(rnorm(n*par.pi[2], 1, sd=par.err[2]), rnorm(n*par.pi[2], 0, sd=par.err[2]))
  X3 = cbind(rnorm(n*par.pi[3], 0, sd=par.err[3]), rnorm(n*par.pi[3], 1, sd=par.err[3]))
  X4 = cbind(rnorm(n*par.pi[4], 1, sd=par.err[4]), rnorm(n*par.pi[4], 1, sd=par.err[4]))
  X = rbind(X1,X2,X3,X4)
  return(X)
}




