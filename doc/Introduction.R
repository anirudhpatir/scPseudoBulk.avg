## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=F, warning=F----------------------------------------------
library(Matrix)
library(scPseudoBulk.avg)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(Seurat)
library(SeuratData)

## ---- message=F, warning=F----------------------------------------------------
InstallData("pbmc3k")
pbmc3k = LoadData("pbmc3k")

## ---- message=F, warning=F, fig.height=5, fig.width=7-------------------------
pbmc3k = SCTransform(pbmc3k, verbose = F)
pbmc3k = RunPCA(pbmc3k, verbose = F)
npca = elbow_knee(pbmc3k@reductions$pca@stdev)
pbmc3k = FindNeighbors(pbmc3k, dims = 1:npca, verbose = F)
pbmc3k = FindClusters(pbmc3k, verbose = F)
ElbowPlot(pbmc3k, ndims = 50) + geom_vline(xintercept = npca)
pbmc3k = RunUMAP(pbmc3k, dims = 1:npca, verbose = F)
DimPlot(pbmc3k, label = T)

## -----------------------------------------------------------------------------
scbulk = scPseudoBulk.avg(pbmc3k@assays$SCT@data, clusters = pbmc3k@active.ident)
head(scbulk)

## ---- fig.height=5, fig.width=7-----------------------------------------------
pheatmap(scbulk, scale = "row", clustering_method = "ward.D2")

