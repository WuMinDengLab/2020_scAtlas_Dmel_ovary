library(RANN)
library(Seurat)
library(ggpubr)
library(dplyr)
library(R.utils)
library(clustree)
library(monocle)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(corrplot)
library(org.Dm.eg.db)
library(topGO)


data_dir<-"U:/Deep/RNA_Seq/scRNA_Seq/Data/FC/round2/Ctrl_FC/"
ctrl.data<-Read10X(data.dir=data_dir)
FC.control <- CreateSeuratObject(raw.data=ctrl.data, min.cells=3, min.genes=150, project="WT Repl1")

#Initial QC sanity checks
ggplot(FC.control@meta.data, aes(nGene)) + geom_histogram(bins=100, color = "black", fill = "white") + xlab("Sum of expression") + ggtitle("Total expression histogram before normalization")
ggplot(FC.control@meta.data, aes(nGene)) + geom_density(alpha = 0.2, color = "red", fill = "red") + xlab("Sum of expression") + ggtitle("Total expression density before normalization")
ggplot(FC.control@meta.data, aes(nUMI)) + geom_histogram(bins=100, color = "black", fill = "white") + xlab("Total nUMI") + ggtitle("Total UMIs histogram before normalization")
ggplot(FC.control@meta.data, aes(nUMI)) + geom_density(alpha = 0.2, color = "red", fill = "red") + xlab("Total nUMI") + ggtitle("Total UMIs density before normalization")
ggplot(FC.control@meta.data, aes(nGene, nUMI)) + geom_point() + geom_smooth(method=lm) + stat_cor(method = "pearson", label.x = 4000, label.y = 3000) + xlab("nGene") + ylab("nUMI")

mito.genes <- grep(pattern = "-m", x = rownames(x = FC.control@data), value = TRUE)
percent.mito <- Matrix::colSums(FC.control@raw.data[mito.genes, ])/Matrix::colSums(FC.control@raw.data)
FC.control <- AddMetaData(object = FC.control, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = FC.control, features.plot = c("nGene", "nUMI", "percent.mito"), point.size.use = 0.8, nCol = 3)
table(FC.control@meta.data$orig.ident)
FC.control <- FilterCells(object = FC.control, subset.names = c("nGene", "nUMI", "percent.mito"), low.thresholds = c(775, -Inf, -Inf), high.thresholds = c(2200, 18000, 0.005))
FC.control <- NormalizeData(object = FC.control, normalization.method = "LogNormalize", scale.factor = 10000)
ggplot(FC.control@meta.data, aes(nGene)) + geom_histogram(bins=100, color = "black", fill = "white") + xlab("Sum of expression") + ggtitle("Total expression histogram after normalization")
ggplot(FC.control@meta.data, aes(nGene)) + geom_density(alpha = 0.2, color = "red", fill = "red") + xlab("Sum of expression") + ggtitle("Total expression density after normalization")
ggplot(FC.control@meta.data, aes(nUMI)) + geom_histogram(bins=100, color = "black", fill = "white") + xlab("Total nUMI") + ggtitle("Total UMIs histogram after normalization") + xlim(0,20000)
ggplot(FC.control@meta.data, aes(nUMI)) + geom_density(alpha = 0.2, color = "red", fill = "red") + xlab("Total nUMI") + ggtitle("Total UMIs density after normalization") + xlim(0,20000)
table(FC.control@meta.data$orig.ident)
VlnPlot(object = FC.control, ident.include = FC.control@meta.data$orig.ident, features.plot = c("nGene", "nUMI", "percent.mito"), point.size.use = 0.8, nCol = 3)
ggplot(FC.control@meta.data, aes(nGene, nUMI)) + geom_point() + geom_smooth(method=lm) + stat_cor(method = "pearson", label.x = 1500, label.y = 3000) + xlab("nGene") + ylab("nUMI")
length(x = FC.control@var.genes)
table(FC.control@meta.data$orig.ident)

#Post Normalization QC sanity checks
x <- summary(FC.control@meta.data$nGene)
x <- FC.control@meta.data$nGene
x_1 <- (x-mean(x))/sd(x)
summary(x_1)
df <- data.frame(x=x_1)
ggplot(df, aes(x)) + geom_histogram(bins=100, color = "black", fill = "white") + ggtitle("Distribution of total expression after normalization") + geom_vline(xintercept = c(-2,2), linetype="dotted", color="red")

#CellCycle Scoring
cc.genes <- readLines(con = "U:/Deep/RNA_Seq/scRNA_Seq/cell_cycle_genes.txt")
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
DimHeatmap(object = FC.control, reduction.type = "pca", cells.use = 500, dim.use = 61:75, do.balanced = TRUE)

PCElbowPlot(FC.control, num.pc = 100)

use.pcs = c(1:75)
PCElbowPlot(FC.control, num.pc = 100) + geom_vline(xintercept = c(0,75), linetype="dashed", color="red")
FC.control <- FindClusters(object = FC.control, reduction.type = "pca", dims.use = use.pcs, resolution = 4.5, print.output = FALSE, save.SNN = TRUE, algorithm = 3)
PrintFindClustersParams(object = FC.control)
sapply(grep("^res",colnames(FC.control@meta.data),value = TRUE), function(x) length(unique(FC.control@meta.data[,x])))
#clustree(FC.control, layout = "sugiyama", use_core_edges = FALSE)

FC.control <- SetAllIdent(FC.control, id = "res.4.5")
table(FC.control@ident,FC.control@meta.data$orig.ident)
#FC.control <- RunTSNE(object = FC.control, reduction.use = "pca", dims.use = use.pcs, do.fast = TRUE)
FC.control <- RunUMAP(object = FC.control, reduction.use = "pca", dims.use = use.pcs, n_neighbors = 20, min_dist = 0.35)
DimPlot(object = FC.control, reduction = "umap", pt.size = 2, do.label = TRUE)

#Merging Similar Clusters
FC.control <- RenameIdent(object = FC.control, old.ident.name = c('11'), new.ident.name = '0')
FC.control <- RenameIdent(object = FC.control, old.ident.name = c('20'), new.ident.name = '0')
FC.control <- RenameIdent(object = FC.control, old.ident.name = c('22'), new.ident.name = '0')
FC.control <- RenameIdent(object = FC.control, old.ident.name = c('17'), new.ident.name = '0')
FC.control <- RenameIdent(object = FC.control, old.ident.name = c('13'), new.ident.name = '0')
FC.control <- RenameIdent(object = FC.control, old.ident.name = c('25'), new.ident.name = '0')

FC.control <- RenameIdent(object = FC.control, old.ident.name = c('14'), new.ident.name = '9')
FC.control <- RenameIdent(object = FC.control, old.ident.name = c('34'), new.ident.name = '9')
FC.control <- RenameIdent(object = FC.control, old.ident.name = c('15'), new.ident.name = '9')

FC.control <- RenameIdent(object = FC.control, old.ident.name = c('26'), new.ident.name = '10')

FC.control <- RenameIdent(object = FC.control, old.ident.name = c('28'), new.ident.name = '7')

DimPlot(object = FC.control, reduction = "umap", pt.size = 1, do.label = TRUE)

#MANUALLY SELECT GERMLINE CLUSTER (8,23) FOR FILTERING AND IDENTIFICATION:
select.cells <- DimPlot(object = FC.control, reduction = "umap", pt.size = 1, do.label = TRUE, do.identify = TRUE)
head(select.cells)
FC.control <- SetIdent(object = FC.control, cells.use = select.cells, ident.use = "8")
FC.control <- SubsetData(FC.control, ident.remove = c(23))

clueless <- StashIdent(object = clueless, save.name = "CellType")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = clueless, ident = "0"), ident.use = "Contamination")

DimPlot(object = FC.control, reduction = "umap", pt.size = 1, do.label = FALSE)


top5 <- markers_all_single %>% group_by(cluster) %>% top_n(5, avg_logFC)
dim(top5)
DoHeatmap(object = FC.control, genes.use = top5$gene, slim.col.label = TRUE, remove.key = TRUE)
tiff("Plots/Hmap.tiff", height = 8000, width = 10000, res = 500)
DoHeatmap(object = FC.control, genes.use = top5$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()


DoHeatmap(object = clueless, genes.use = top15$gene, slim.col.label = TRUE, remove.key = FALSE, col.low = "black", col.mid = "red4", col.high = "yellow")

FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Hemocytes"), ident.use = "Hemocytes")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Fat Body Cells"), ident.use = "Adipocytes")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Oviduct 2"), ident.use = "Oviduct Cells 2")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Oviduct 1"), ident.use = "Oviduct Cells 1")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Muscle Cells 2"), ident.use = "Muscle Sheath Cells 2")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Muscle Cells 1"), ident.use = "Muscle Sheath Cells 1")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Floor Cells"), ident.use = "Dorsal Appendage Cells 2")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Roof Cells"), ident.use = "Dorsal Appendage Cells 1")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Migrating Border FC"), ident.use = "Anterior Follicle Cells 2")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Nurse Cell Associating FC 3"), ident.use = "Anterior Follicle Cells 4")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Nurse Cell Associating FC 2"), ident.use = "Anterior Follicle Cells 3")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Nurse Cell Associating FC 1"), ident.use = "Anterior Follicle Cells 1")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Corpus luteum 2"), ident.use = "Anterior Follicle Cells 5")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Corpus luteum 1"), ident.use = "Corpus luteum Cells")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Oocyte Associating FC 6"), ident.use = "Post-Mitotic Main Body Follicle Cells 8")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Oocyte Associating FC 5"), ident.use = "Post-Mitotic Main Body Follicle Cells 7")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Oocyte Associating FC 4"), ident.use = "Post-Mitotic Main Body Follicle Cells 6")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Oocyte Associating FC 3"), ident.use = "Post-Mitotic Main Body Follicle Cells 5")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Oocyte Associating FC 2"), ident.use = "Post-Mitotic Main Body Follicle Cells 4")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Oocyte Associating FC 1"), ident.use = "Post-Mitotic Main Body Follicle Cells 3")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Endocycling FC 2"), ident.use = "Post-Mitotic Main Body Follicle Cells 2")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Endocycling FC 1"), ident.use = "Post-Mitotic Main Body Follicle Cells 1")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Transitioning FC"), ident.use = "Mitotic-Endocycle Transitioning Follicle Cells")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Mitotic FC"), ident.use = "Mitotic Follicle Cells")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Stalk Cells"), ident.use = "Stalk Cells")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Polar Cells"), ident.use = "Polar Cells")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Somatic Cells of the Germline Niche"), ident.use = "Somatic Cells of the Germarium")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Late Germline Cells"), ident.use = "Late Germline Cells")
FC.control <- SetIdent(object = FC.control, cells.use = WhichCells(object = FC.control, ident = "Early Germline Cells"), ident.use = "Early Germline Cells")

FC.avg=AverageExpression(FC.control, use.scale = TRUE, show.progress = TRUE)
plot(FC.avg[,"Muscle Cells 1"],FC.avg[,"Muscle Cells 2"],pch=20,cex=1,xlab="Muscle Cells 1",ylab="Muscle Cells 2",main=paste0("Muscle Cells 1 vs 2 (cor= ", formatC(cor(FC.avg[,"Muscle Cells 1"],FC.avg[,"Muscle Cells 2"]), 2, format="f"),")")) + abline(lm(FC.avg[,"Muscle Cells 2"]~FC.avg[,"Muscle Cells 1"]), col="red")


FC.subset <- SubsetData(FC.control, ident.remove = c("Early Germline Cells", "Late Germline Cells", "Fat Body Cells", "Hemocytes", "Oviduct 1", "Oviduct 2"))

FC.subset <- NormalizeData(object = FC.subset, normalization.method = "LogNormalize", scale.factor = 10000)
FC.avg=AverageExpression(FC.subset, use.scale = TRUE, show.progress = TRUE)
M <- cor(FC.avg)
head(M)
head(round(M,2))
res1 <- cor.mtest(M, conf.level = .95)
tiff("Plots/Correlation.tiff", height = 10000, width = 8000, res = 800)
corrplot(M, method = "color", p.mat = res1$p, pch.cex = 1.2, addrect = 5, insig = "label_sig", pch.col = "black", tl.col = "black", order="hclust", col=brewer.pal(n=9, name="YlOrRd"), outline = "white")
dev.off()


tiff("Plots/DotPlot.tiff", height = 5800, width = 8600, res = 800)
features.plot <- c("Hml", "Ilp6", "abd-A", "Mp20", "Ance", "Mmp2", "Glut4EF", "peb", "Past1", "slbo", "Cad99C", "mirr", "Diap1", "Ilp8", "Femcoat", "yellow-g", "Cad74A", "Cad86C", "Fcp3C", "ttk", "dec-1", "Vml", "psd", "ct", "Atf3", "upd1", "Wnt4", "orb", "vas", "bam")
DotPlot(object = FC.control, genes.plot = features.plot, plot.legend = TRUE, x.lab.rot = T, dot.scale = 6, cols.use = c("gold", "red"))
dev.off()


FC.control <- SetIdent(FC.control, cells.use = WhichCells(FC.control,c("1. Early Germline Cells", "2. Late Germline Cells", "3. Somatic Cells in Germarium Region 1", "4. Polar Follicle Cells", "5. Interfollicle (Stalk) Cells", "6. Pre-Mitotic Follicle Cells", "7. Mitotic Follicle Cells", "10. Post-Mitotic Main Body Follicle Cells 1", "11. Post-Mitotic Main Body Follicle Cells 2", "12. Post-Mitotic Main Body Follicle Cells 3", "13. Post-Mitotic Main Body Follicle Cells 4", "14. Post-Mitotic Main Body Follicle Cells 5", "15. Post-Mitotic Main Body Follicle Cells 6", "Stage 10B Main Body Follicle Cells", "Stage 12 FC (yellow-g)", "Dorsal Appendage forming FC (Cad86C)", "Stage 14 FC (Femcoat)", "Stage 14 FC (Men/Ilp8)", "20. Anterior Follicle Cells 2", "21. Anterior Follicle Cells 3", "22. Anterior Follicle Cells 4", "23. Anterior Follicle Cells 5", "24. Main Body-derived Corpus luteum Cells", "Mmp2+ Terminal CL Cells", "Terminal CL Cells (Ance/Mmp2)", "27. Muscle Sheath Cells", "28. Oviduct Cells", "29. Adipocytes", "30. Hemocytes")), ident.use = "Else")
FC.control <- SetIdent(FC.control, cells.use = WhichCells(FC.control, c("8. Post-Mitotic Follicle Cells 1", "9. Post-Mitotic Follicle Cells 2", "19. Anterior Follicle Cells 1")), ident.use = "Subset")
DimPlot(object = FC.control, reduction = "umap", pt.size = 1, do.label = FALSE, cols.use = c("blue", "firebrick2"))
tiff("Plots/umap_Figure_5.tiff", height = 6000, width = 6000, res = 800)
DimPlot(object = FC.control, reduction = "umap", pt.size = 1.4, cols.use = c("lightgrey", "blue"))
dev.off()


#Germline Trajectory
clueless_monocle <- SubsetData(FC.control, ident.use = c("Germline Stem Cells", "Cystoblasts", "clueless"))
clueless_cds <- importCDS(clueless_monocle, import_all = TRUE)
clueless_cds <- estimateSizeFactors(clueless_cds)
clueless_cds <- estimateDispersions(clueless_cds)
clueless_cds <- detectGenes(clueless_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(clueless_cds), num_cells_expressed >= 20))
#plot_pc_variance_explained(clueless_cds, return_all = F)
clueless_cds<- reduceDimension(clueless_cds, max_components = 2, norm_method = 'log', num_dim = 7, reduction_method = 'tSNE', verbose = T)
clueless_cds <- clusterCells(clueless_cds, verbose = F)

plot_rho_delta(clueless_cds, rho_threshold = 2, delta_threshold = 4) + xlab("Local density of cluster (rho)") + ylab("Cell-Cell Distance (delta)")
clueless_cds <- clusterCells(clueless_cds, rho_threshold = 10, delta_threshold = 5, skip_rho_sigma = T, verbose = F)
plot_rho_delta(clueless_cds, rho_threshold = 2, delta_threshold = 4) + xlab("Local density of cluster (rho)") + ylab("Cell-Cell Distance (delta)") + annotate("rect", xmin = 10, xmax = 16, ymin = 5, ymax = 40, alpha = .2, color = "red", fill = "red")
table(pData(clueless_cds)$Cluster)
clustering_DEG_genes <- differentialGeneTest(clueless_cds[expressed_genes,], fullModelFormulaStr = '~Cluster', cores = 8)

ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:200]
clueless_cds <- setOrderingFilter(clueless_cds, ordering_genes = ordering_genes)
plot_ordering_genes(clueless_cds)
clueless_cds <- reduceDimension(clueless_cds, method = 'DDRTree')
clueless_cds <- orderCells(clueless_cds)

plot_cell_trajectory(clueless_cds, color_by = "Pseudotime")
#To reverse the pseudotime START/END POINT
clueless_cds <- orderCells(clueless_cds, reverse = TRUE)

#Plot Clusters on tSNE space
plot_cell_clusters(clueless_cds, cell_size = 3) +
    theme(legend.text=element_text(size=10)) + #set the size of the text
    theme(legend.position="right") #put the color legend on the right
#Plot Clusters on a trajectory
plot_cell_trajectory(clueless_cds, color_by = "Cluster") +
    theme(legend.text=element_text(size=10)) + #set the size of the text
    theme(legend.position="right") #put the color legend on the right
#Identify different States on a trajectory
plot_cell_trajectory(clueless_cds, color_by = "State") +
    theme(legend.text=element_text(size=10)) + #set the size of the text
    theme(legend.position="right") #put the color legend on the right

plot_cell_clusters(clueless_cds, cell_size = 3, color_by = "Pseudotime") +
    theme(legend.text=element_text(size=10)) + #set the size of the text
    theme(legend.position="right") + #put the color legend on the right
    scale_color_viridis(option = "B")





clueless_monocle <- SubsetData(FC.control, ident.use = c("Escort Cells", "Polar Cells", "Stalk Cells"))
clueless_cds <- importCDS(clueless_monocle, import_all = TRUE)
clueless_cds <- estimateSizeFactors(clueless_cds)
clueless_cds <- estimateDispersions(clueless_cds)
clueless_cds <- detectGenes(clueless_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(clueless_cds), num_cells_expressed >= 50))
#plot_pc_variance_explained(clueless_cds, return_all = F)
clueless_cds<- reduceDimension(clueless_cds, max_components = 2, norm_method = 'log', num_dim = 16, reduction_method = 'tSNE', verbose = T)
clueless_cds <- clusterCells(clueless_cds, verbose = F)

plot_rho_delta(clueless_cds, rho_threshold = 2, delta_threshold = 4) + xlab("Local density of cluster (rho)") + ylab("Cell-Cell Distance (delta)")
clueless_cds <- clusterCells(clueless_cds, num_clusters = 7, verbose = F)
table(pData(clueless_cds)$Cluster)

clustering_DEG_genes <- differentialGeneTest(clueless_cds[expressed_genes,], fullModelFormulaStr = '~CellType', cores = 8)

ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
clueless_cds <- setOrderingFilter(clueless_cds, ordering_genes = ordering_genes)
plot_ordering_genes(clueless_cds)
clueless_cds <- reduceDimension(clueless_cds, method = 'DDRTree')
clueless_cds <- orderCells(clueless_cds)
plot_cell_clusters(clueless_cds, color_by = "CellType", cell_size = 3)
plot_cell_trajectory(clueless_cds, color_by = "CellType", cell_size = 2)

plot_cell_trajectory(clueless_cds, color_by = "Pseudotime", markers = c("mirr", "Atf3", "LamC", "ct"), use_color_gradient = T)
......................................................................................................................................................

#EARLY SOMATIC
expressed_genes <- row.names(subset(fData(clueless_cds), num_cells_expressed >= 20))
clueless_cds <- clusterCells(clueless_cds, num_clusters = 8, verbose = F)
table(pData(clueless_cds)$Cluster)
clustering_DEG_genes <- differentialGeneTest(clueless_cds[expressed_genes,], fullModelFormulaStr = '~Cluster', cores = 8)

ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:895]
clueless_cds <- setOrderingFilter(clueless_cds, ordering_genes = ordering_genes)
clueless_cds <- reduceDimension(clueless_cds, method = 'DDRTree')
clueless_cds <- orderCells(clueless_cds)
plot_cell_clusters(clueless_cds, color_by = "Cluster", cell_size = 3)
plot_cell_trajectory(clueless_cds, color_by = "Cluster", cell_size = 2)

my_vector <- rep('no', nrow(pData(clueless_cds)))
my_vector[pData(clueless_cds)$Cluster == 1] <- rep('Polar Cells', sum(pData(clueless_cds)$Cluster == 1))
my_vector[pData(clueless_cds)$Cluster == 2] <- rep('Somatic Stem Cells', sum(pData(clueless_cds)$Cluster == 2))
my_vector[pData(clueless_cds)$Cluster == 3] <- rep('Mitotic Follicle Cells', sum(pData(clueless_cds)$Cluster == 3))
my_vector[pData(clueless_cds)$Cluster == 4] <- rep('Stalk Cells', sum(pData(clueless_cds)$Cluster == 4))
my_vector[pData(clueless_cds)$Cluster == 5] <- rep('Mitotic Follicle Cells', sum(pData(clueless_cds)$Cluster == 5))
my_vector[pData(clueless_cds)$Cluster == 6] <- rep('Mitotic Follicle Cells', sum(pData(clueless_cds)$Cluster == 6))
my_vector[pData(clueless_cds)$Cluster == 7] <- rep('Mitotic Follicle Cells', sum(pData(clueless_cds)$Cluster == 7))
my_vector[pData(clueless_cds)$Cluster == 8] <- rep('Mitotic Follicle Cells', sum(pData(clueless_cds)$Cluster == 8))
my_vector[pData(clueless_cds)$Cluster == 9] <- rep('pre-Follicle Cells', sum(pData(clueless_cds)$Cluster == 9))
my_vector[pData(clueless_cds)$Cluster == 10] <- rep('Pre-Follicle Cells', sum(pData(clueless_cds)$Cluster == 10))
pData(clueless_cds)$CellType <- my_vector

my_genes <- row.names(subset(fData(clueless_cds),
                             gene_short_name %in% c("Alh", "how", "E(spl)mbeta-HLH", "dpn", "CtBP", "lilli", "Kdm4B", "Glut4EF", "cas", "zfh1", "Atf3", "Trf2", "aop", "sd", "Eip75B", "nvy", "Nf-YA", "Sox14", "E(spl)m7-HHLH", "nclb", "E(y)3", "dj-1beta", "Jra", "eya", "HmgD", "h", "Nacalpha")))
tiff("TF_somatic_Hmap.tiff", height = 2000, width = 3000, res = 800)
plot_genes_branched_heatmap(clueless_cds[my_genes,], norm_method = "log", branch_point = 1, num_clusters = 1, cores = 8, use_gene_short_name = TRUE, show_rownames = TRUE)
dev.off()

my_genes <- row.names(subset(fData(clueless_cds),
                             gene_short_name %in% c("wdb", "Pka-R2", "dpn", "His2Av", "E(spl)m7-HLH", "Gbeta13F", "Rac2", "mav", "nmo", "spz", "18w", "dl", "vih", "eff", "sina", "shg", "fz", "fzy", "stg", "Pura", "Cad99C", "Shc", "14-3-3zeta", "hh", "Wnt4", "Wnt6", "FK506-bp1", "FK506-bp2", "dUTPase", "His2Av", "sd", "N")))
tiff("SP_somatic_Hmap.tiff", height = 2500, width = 3000, res = 800)
plot_genes_branched_heatmap(clueless_cds[my_genes,], norm_method = "log", branch_point = 1, num_clusters = 1, cores = 8, use_gene_short_name = TRUE, show_rownames = TRUE)
dev.off()

my_genes <- row.names(subset(fData(clueless_cds),
                             gene_short_name %in% c("sog", "Med", "Fur1", "Rac2", "rl", "Cbl", "nmo", "Dsor1", "Egfr", "wbd", "sl", "14-3-3zeta", "CkIalpha", "Pka-R2", "E(spl)mbeta-HLH", "dpn", "E(spl)m3-HLH", "fng", "dl", "ps", "E(spl)m7-HLH", "Gbeta13F", "sgg", "Tak1", "Spn27A", "spz", "fz2", "CtBP", "Pi3K21B")))
tiff("SP_somatic_Hmap.tiff", height = 2500, width = 3000, res = 800)
plot_genes_branched_heatmap(clueless_cds[my_genes,], norm_method = "log", branch_point = 1, num_clusters = 1, cores = 8, use_gene_short_name = TRUE, show_rownames = TRUE)
dev.off()



tiff("Plots/early_somatic_TFState_Cluster.tiff", height = 450, width = 550, res = 100)
plot_cell_clusters(clueless_cds, color_by = "State", cell_size = 3) +
    theme(legend.text=element_text(size=10)) + #set the size of the text
    theme(legend.position="right") #put the color legend on the right
dev.off()
tiff("Plots/early_somatic_Pseudotime_Cluster.tiff", height = 450, width = 610, res = 100)
plot_cell_clusters(clueless_cds, color_by = "Pseudotime", cell_size = 3) +
    theme(legend.text=element_text(size=10)) + #set the size of the text
    theme(legend.position="right") #put the color legend on the right
dev.off()
tiff("Plots/early_somatic_CellType_Cluster.tiff", height = 450, width = 650, res = 100)
plot_cell_clusters(clueless_cds, color_by = "CellType", cell_size = 3) +
    theme(legend.text=element_text(size=10)) + #set the size of the text
    theme(legend.position="right") #put the color legend on the right
dev.off()

genes <- row.names(subset(fData(clueless_cds),
                               gene_short_name %in% c("bbg", "CG14207", "Drat", "emc", "trio", "stck", "CG1620", "stai", "hdc", "cnc", "eya", "cas", "Nrx-IV", "sli", "HmgD", "D1")))
tiff("Plots/boopboop.tiff", height = 4000, width = 6000, res = 800)
plot_genes_branched_pseudotime(clueless_cds[genes,],
                               branch_point = 1,
                               color_by = "Pseudotime", cell_size = 1,
                               ncol = 4, trend_formula = "~ sm.ns(Pseudotime, df=2) * Branch")
dev.off()

......................................................................................................................................................
#POLAR AND STALK
ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][c(2:19,23:35)]
my_cds_subset <- setOrderingFilter(my_cds_subset, ordering_genes = ordering_genes)
my_cds_subset <- reduceDimension(my_cds_subset, method = 'DDRTree')
my_cds_subset <- orderCells(my_cds_subset)
plot_cell_trajectory(my_cds_subset, color_by = "Cluster", cell_size = 2)

my_genes <- row.names(subset(fData(my_cds_subset),
                             gene_short_name %in% c("nmo", "Wdr62", "sona", "CG45263", "jbug", "CG14223", "sli", "mbl", "kek1", "bbg", "mfas", "CG43658", "upd1", "Hk", "CG6357", "E(spl)m7-HLH", "GILT1", "CG2157", "CR43930", "Jupiter", "Fas3", "mamo", "RapGAP1", "drongo", "rgn", "Fas2", "CG18208", "chic", "CG9336", "CG1572", "CG1648", "GILT2", "CG15211", "GILT3", "GILT1", "CG10311", "CG13377", "ringer", "CG13403", "CG42797", "bbg", "E(spl)m2-BFM", "brat", "Fas2", "CG43858", "Alh", "sona", "CG1648", "LamC")))
tiff("Polar_Stalk_Hmap.tiff", height = 5000, width = 3000, res = 800)
plot_genes_branched_heatmap(my_cds_subset[my_genes,], norm_method = "log", branch_point = 1, num_clusters = 2, cores = 8, use_gene_short_name = TRUE, show_rownames = TRUE)
dev.off()
......................................................................................................................................................
#GERMARIUM
ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][c(2:41,43:73,75:100)]
stem_cds <- setOrderingFilter(stem_cds, ordering_genes = ordering_genes)
plot_ordering_genes(stem_cds)
stem_cds <- reduceDimension(stem_cds, method = 'DDRTree')
stem_cds <- orderCells(stem_cds)
plot_cell_clusters(stem_cds, color_by = "Cluster", cell_size = 3)
plot_cell_trajectory(stem_cds, color_by = "Cluster", cell_size = 2)

marker_genes <- row.names(subset(fData(clueless_cds), gene_short_name %in% c("Nprl2","zf30C","rgr","Six4","D19A","Smox","ear","D12","ftz-f1","NFAT","pzg","Stat92E","Ssl1","Mms19","nej","row","Tfb1","dpn","CG9890","Max","gro","cg","so","ken","CtBP","Sin3A","e(y)3","l(3)mbt","ush","cas","Hr96","nclb","l(1)10Bb","Tif-IA","CG2199","Dp","H15","Aatf","sob","sr","Alh","mirr","Rpn11","Tudor-SN","Hr39","HmgZ","CG5382","Hnf4","slbo","Hr4","Eip74EF","l(3)neo38","Trs20","vri","tj", "Sox14", "spt4")))

diff_test_res <- differentialGeneTest(clueless_cds[marker_genes,],fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 8)
plot_pseudotime_heatmap(clueless_cds[row.names(subset(diff_test_res, qval < 0.05 & num_cells_expressed > 50)),], num_clusters = 3,
                                                       cores = 8,
                                                       show_rownames = T)


my_genes <- row.names(subset(fData(clueless_cds),
                             gene_short_name %in% c("Vps4",
                                                    "betaTub56D",
                                                    "Arpc4",
                                                    "ctp",
                                                    "srp",
                                                    "twi")))
tiff("Plots/EndoMorpho_Hmap.tiff", height = 4000, width = 3000, res = 800)
plot_genes_branched_heatmap(clueless_cds[my_genes,], norm_method = "log", branch_point = 1, num_clusters = 2, cores = 8, use_gene_short_name = TRUE, show_rownames = TRUE)
dev.off()


ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:775]
clueless_cds <- setOrderingFilter(clueless_cds, ordering_genes = ordering_genes)
clueless_cds <- reduceDimension(clueless_cds, method = 'DDRTree')
clueless_cds <- orderCells(clueless_cds)
plot_cell_trajectory(clueless_cds, color_by = "Cluster", cell_size = 2)

BEAM_res <- BEAM(clueless_cds, branch_point = 1, cores = 8)
					BEAM_res <- BEAM_res[-grep("^Yp",BEAM_res$gene_short_name),]
					BEAM_res <- BEAM_res[-grep("^Cp",BEAM_res$gene_short_name),]
					BEAM_res <- BEAM_res[-grep("^Vm",BEAM_res$gene_short_name),]
					BEAM_res <- BEAM_res[-grep("^CR",BEAM_res$gene_short_name),]
					BEAM_res <- BEAM_res[-grep("^Obp",BEAM_res$gene_short_name),]
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval", "num_cells_expressed")]
tiff("Plots/CL_Hmap.tiff", height = 45000, width = 2000, res = 400)
clueless_branched_heatmap <- plot_genes_branched_heatmap(clueless_cds[row.names(subset(BEAM_res, qval < 0.05 & num_cells_expressed > 50)),], norm_method = "log", branch_point = 1, num_clusters = 5, cores = 8, use_gene_short_name = TRUE, show_rownames = TRUE, return_heatmap = TRUE)
dev.off()

