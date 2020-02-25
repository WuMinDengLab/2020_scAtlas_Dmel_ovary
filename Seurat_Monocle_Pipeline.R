library(Seurat)
library(ggpubr)
library(dplyr)
library(clustree)
library(monocle)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(corrplot)

## SEURAT 2.3.4 PIPELINE ##

data_dir<-"~/cellranger_output/"
ctrl.data<-Read10X(data.dir=data_dir)
FC.control <- CreateSeuratObject(raw.data=ctrl.data, min.cells=3, min.genes=775, project="WT Repl1")

mito.genes <- grep(pattern = "-m", x = rownames(x = FC.control@data), value = TRUE)
percent.mito <- Matrix::colSums(FC.control@raw.data[mito.genes, ])/Matrix::colSums(FC.control@raw.data)
FC.control <- AddMetaData(object = FC.control, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = FC.control, features.plot = c("nGene", "nUMI", "percent.mito"), point.size.use = 0.8, nCol = 3)
FC.control <- FilterCells(object = FC.control, subset.names = c("nGene", "nUMI", "percent.mito"), low.thresholds = c(775, -Inf, -Inf), high.thresholds = c(2200, 18000, 0.01))
FC.control <- NormalizeData(object = FC.control, normalization.method = "LogNormalize", scale.factor = 10000)
table(FC.control@meta.data$orig.ident)

#CellCycle Scoring
cc.genes <- readLines(con = "~/cell_cycle_genes.txt") ##PROVIDED IN THE SUPPLEMENTARY
g2m.genes <- cc.genes[1:68]
s.genes <- cc.genes[69:124]
FC.control <- CellCycleScoring(object = FC.control, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)
head(x = FC.control@meta.data)
FC.control@meta.data$CC.Difference <- FC.control@meta.data$S.Score - FC.control@meta.data$G2M.Score
FC.control <- ScaleData(object = FC.control, vars.to.regress = c("nUMI", "percent.mito", "CC.Difference"), display.progress = TRUE)
FC.control <- FindVariableGenes(object = FC.control, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff= 0.5)
length(x = FC.control@var.genes)

#PCA
FC.control <- RunPCA(object = FC.control, pc.genes = FC.control@var.genes, do.print = TRUE, pcs.compute = 100, pcs.print = 1:5, genes.print = 10)
PCElbowPlot(FC.control, num.pc = 100)

#Choose number of PCs
use.pcs = c(1:30)
PCElbowPlot(FC.control, num.pc = 100) + geom_vline(xintercept = c(0,75), linetype="dashed", color="red")
FC.control <- FindClusters(object = FC.control, reduction.type = "pca", dims.use = use.pcs, resolution = 4.5, print.output = FALSE, save.SNN = TRUE, algorithm = 3)
PrintFindClustersParams(object = FC.control)
sapply(grep("^res",colnames(FC.control@meta.data),value = TRUE), function(x) length(unique(FC.control@meta.data[,x])))

#Compare different resolutions
clustree(FC.control, layout = "sugiyama", use_core_edges = FALSE)

#Choose resolution
FC.control <- SetAllIdent(FC.control, id = "res.4.5")
table(FC.control@ident,FC.control@meta.data$orig.ident)

#Run UMAP and plot
FC.control <- RunUMAP(object = FC.control, reduction.use = "pca", dims.use = use.pcs, n_neighbors = 20, min_dist = 0.35)
DimPlot(object = FC.control, reduction = "umap", pt.size = 1.2, do.label = TRUE)



## MONOCLE v2 PIPELINE ##

seurat2monocle <- SubsetData(FC.control, ident.use = c("###")) ### Choose Seurat Clusters to Subset
trajectory_cds <- importCDS(seurat2monocle, import_all = TRUE)
trajectory_cds <- estimateSizeFactors(trajectory_cds)
trajectory_cds <- estimateDispersions(trajectory_cds)
trajectory_cds <- detectGenes(trajectory_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(trajectory_cds), num_cells_expressed >= 20))
plot_pc_variance_explained(trajectory_cds, return_all = F)
trajectory_cds<- reduceDimension(trajectory_cds, max_components = 2, norm_method = 'log', num_dim = "###", reduction_method = 'tSNE', verbose = T) ### Choose appropriate dimensions
trajectory_cds <- clusterCells(trajectory_cds, verbose = F)

plot_rho_delta(trajectory_cds, rho_threshold = "###", delta_threshold = "###") + xlab("Local density of cluster (rho)") + ylab("Cell-Cell Distance (delta)") ### Choose appropriate thresholds
trajectory_cds <- clusterCells(trajectory_cds, rho_threshold = 10, delta_threshold = 5, skip_rho_sigma = T, verbose = F)
table(pData(trajectory_cds)$Cluster)
clustering_DEG_genes <- differentialGeneTest(trajectory_cds[expressed_genes,], fullModelFormulaStr = '~Cluster', cores = 8) ### Can also be clustered by ~CellType, instead of ~Cluster

ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)]["###"] ### Evaluate gene list and remove contaminants
trajectory_cds <- setOrderingFilter(trajectory_cds, ordering_genes = ordering_genes)
plot_ordering_genes(trajectory_cds)
trajectory_cds <- reduceDimension(trajectory_cds, method = 'DDRTree')
trajectory_cds <- orderCells(trajectory_cds)

plot_cell_trajectory(trajectory_cds, color_by = "Pseudotime")

#Plot Clusters on tSNE space
plot_cell_clusters(trajectory_cds, cell_size = 3) +
    theme(legend.text=element_text(size=10)) + #set the size of the text
    theme(legend.position="right") #put the color legend on the right
#Plot Clusters on a trajectory
plot_cell_trajectory(trajectory_cds, color_by = "Cluster") +
    theme(legend.text=element_text(size=10)) + #set the size of the text
    theme(legend.position="right") #put the color legend on the right
#Identify different States on the trajectory
plot_cell_trajectory(trajectory_cds, color_by = "State") +
    theme(legend.text=element_text(size=10)) + #set the size of the text
    theme(legend.position="right") #put the color legend on the right
#Identify Pseudotemporal direction of the trajectory
plot_cell_clusters(trajectory_cds, cell_size = 3, color_by = "Pseudotime") +
    theme(legend.text=element_text(size=10)) + #set the size of the text
    theme(legend.position="right") + #put the color legend on the right
    scale_color_viridis(option = "B")

