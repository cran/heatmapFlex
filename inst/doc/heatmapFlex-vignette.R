## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(heatmapFlex)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("heatmapFlex")

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("vfey/heatmapFlex")

## ---- out.width='85%', fig.width=6, fig.height=4, fig.align='center'----------
mat <- matrix(rnorm(100), nrow = 10)
heatmap.n2(mat)

## ---- out.width='85%', fig.width=6, fig.height=4, fig.align='center'----------
mat <- matrix(c(rnorm(50, mean = 1), rnorm(50, mean = -1)), nrow = 10)
heatmap.n2(mat, col = "BuWtRd", rowMembers=rep(1:2, each=5),
           colMembers=rep(1:2, each=5),
           labRow=paste0("gene-", 1:10),
           labCol=paste0(c("A", "B"), rep(1:5, 2)), r.cex=0.8,
           dendroheight = lcm(2.2), dendrowidth = lcm(2.4))

## ---- eval=FALSE--------------------------------------------------------------
#  mat <- matrix(c(rnorm(50, mean = 1), rnorm(50, mean = -1)), nrow = 10)
#  dl <- heatmap.n2(mat, col = "BuWtRd", labRow=paste0("gene-", 1:10),
#                   labCol=paste0(c("A", "B"), rep(1:5, 2)),
#                   r.cex=0.8, dendroheight = lcm(2.2), dendrowidth = lcm(2.4))
#  zoom_heatmap(dl)

## ---- out.width='85%', fig.width=6, fig.height=6, fig.align='center'----------
mat <- matrix(c(rnorm(50, mean = 1), rnorm(50, mean = -1)), nrow = 10)
pd <- data.frame(female=c(0,0,1,0,1,1,0,1,0,1), male=c(1,1,0,1,0,0,0,0,1,0),
                 row.names = paste0(c("A", "B"), rep(1:5, 2)),
                 undeclared=c(0,0,0,0,0,0,1,0,0,0))
pd
dl <- heatmap.n2(
  mat,
  col = "BuWtRd",
  rowMembers=rep(1:2, each=5),
  colMembers=rep(1:2, each=5),
  labRow=paste0("gene-", 1:10),
  labCol=paste0(c("A", "B"), rep(1:5, 2)),
  r.cex=0.8,
  dendroheight = lcm(2.2),
  dendrowidth = lcm(2.4),
  sidebars = list(left=data.frame(min=apply(mat, 1, min), max=apply(mat, 1, max)),
                  bottom=data.frame(
                    mean=apply(mat, 2, mean, na.rm=TRUE),
                    treat=factor(rep(c("A", "B"), 5)))),
  factorpalettefn = colorRampPalette(c("lightblue", "limegreen")),
  picketdata = pd)

