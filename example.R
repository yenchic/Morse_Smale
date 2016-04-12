# Author: Yen-Chi Chen, Christopher R. Genovese, Larry Wasserman
# Maintainer: Yen-Chi Chen <ga014528@gmail.com>
# Reference: Statistical Inference using Morse-Smale Complex
# Date: 06/22/2015
source("MSHD.R")
library(mclust) #for GvHD dataset
####
data(GvHD)
D0 = GvHD.control
D1 = GvHD.pos

#### parameters
h = 50  # smoothing
n_g = 21  # grid number
####
D_pull =rbind(D0,D1)

#### grid points
s1 = seq(from=quantile(D_pull[,1], 0.1), to=quantile(D_pull[,1], 0.9), length.out=n_g)
s2 = seq(from=quantile(D_pull[,2], 0.1), to=quantile(D_pull[,2], 0.9), length.out=n_g)
s3 = seq(from=quantile(D_pull[,3], 0.1), to=quantile(D_pull[,3], 0.9), length.out=n_g)
s4 = seq(from=quantile(D_pull[,4], 0.1), to=quantile(D_pull[,4], 0.9), length.out=n_g)
G0 = expand.grid(s1,s2,s3,s4)

#### density estimate
G0_den = kde(D0, G0, h)
G1_den = kde(D1, G0, h)

G_diff = G0_den-G1_den
G_diff_s = G_diff/(max(max(G_diff), abs(min(G_diff))))


#### High dimensional visualization
GvHD_MSHD = MSHD(G0, G_diff_s)


  ### visualization using modes and minima (signatures)
plot.tmp = plot(GvHD_MSHD)
text(plot.tmp$modes.vis, labels = 1:nrow(plot.tmp$modes.vis), cex=3)
text(plot.tmp$minima.vis, labels = 1:nrow(plot.tmp$minima.vis), cex=2)


  ### visualization using cells
plot.tmp = cell.plot(GvHD_MSHD,  r_circle=1, col.cell.filled = "limegreen")
text(plot.tmp$cell.vis, labels = 1:nrow(plot.tmp$cell.vis), cex=2)


  ### visualizing the significant regions
a = quantile(G_diff_s, 0.9)
plot.tmp = cell.plot(GvHD_MSHD, lv.pos=a, lv.neg = -a, r_circle=1, col.cell.filled = "white")
text(plot.tmp$cell.vis, labels = 1:nrow(plot.tmp$cell.vis), cex=2, pos=2)


  ### energy test
E.test = MSEtest(D0,D1,G0)


