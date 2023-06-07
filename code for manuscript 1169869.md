## ## Cellranger(V 4.0)

```
cellranger count --id=run_count_CMV \
   --fastqs=/CMV_fastqs \
   --sample=CMV \
   --transcriptome=/refdata
```

## ## Seurat_SingleR(using R)

#### # packages preparation

```
library(Seurat) 
library(SingleR) 
library(DoubletFinder)
```

#### #Read CellRanger partial matrix files

```
Data <-  Read10X(data.dir="/CMV")
```

#### #Create SeuratObject

```
data <- CreateSeuratObject(counts = Data)
```

#### #Calculate mitochondrial ratio, filter mitochondria, ncount, nfeature

```
data[["percent.mt"]] <- PercentageFeatureSet(data, features=mt_genes)
data <- data@meta.data[which(data@meta.data$percent.mt>=mt), ] %>% nrow()
data <- data@meta.data[which(data@meta.data$nCount_RNA<p_nCount_RNA),] %>% nrow()
data <- data@meta.data[which(data@meta.data$nFeature_RNA<p_nFeature_RNA),] %>% nrow()
```

#### #Remove multi-cells

```
DoubletRate = 0.075 
homotypic.prop <- modelHomotypic(data$seurat_clusters)
nExp_poi <- round(DoubletRate * ncol(data)) 
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
data <- doubletFinder_v3(data, PCs = 1:pcSelect, pN = pN, pK = pK_bcmvn, nExp = nExp_poi.adj)
```

#### #Data normalization

```
data <- NormalizeData(data, normalization.method = method)
```

#### #Analyze highly variable genes

```
data <- FindVariableFeatures(data, selection.method = method, nfeatures = feature)
```

#### #Multi-sample integration

```
Rumen.anchors <- FindIntegrationAnchors(object.list = list, dims = 1:pca) 
Rumen.integrated <- IntegrateData(anchorset = Rumen.anchors, dims = 1:pca)
```

#### #Dimensionality reduction

```
Rumen.integrated <- ScaleData(Rumen.integrated) 
Rumen.integrated <- RunPCA(Rumen.integrated, npcs = pca) 
Rumen.integrated <- RunUMAP(Rumen.integrated, dims = 1:pca) 
Rumen.integrated <- RunTSNE(Rumen.integrated, dims = 1:pca)
```

#### #Generate UMAP plot

```
umap <- DimPlot(Rumen.integrated, reduction = "umap")
```

#### #Generate t-SNE plot

```
tsne <- DimPlot(Rumen.integrated, reduction = "tsne")
```

#### #Clustering

```
Rumen.integrated <- FindNeighbors(Rumen.integrated, dims = 1:pca) 
Rumen.integrated <- FindClusters(Rumen.integrated, resolution = resolutions)
```

#### #Save as rds

```
saveRDS(Rumen.integrated)
```

#### #Differential analysis

```
Rumen.integrated.markers <- FindAllMarkers(Rumen.integrated)
```

#### #Generate heatmap

```
DoHeatmap(Rumen.integrated, features = feature)
```

#### #Generate dot plot

```
DotPlot(Rumen.integrated, features = feature)
```

#### #Read clustering results

```
ctrl <- readRDS(seurat.rds)
```

#### #SingleR annotation

```
test <- as.SingleCellExperiment(ctrl) 
Anno <- SingleR(test = test, ref = BL, method = method)
```

#### #Generate UMAP plot

```
test@meta.data$labels <- Anno$labels
DimPlot(test, reduction = "umap")
```

#### #Generate t-SNE plot

```
DimPlot(test, reduction = "tsne")
```

#### #Generate heatmap

```
DoHeatmap(test, features = feature)
```

#### #Generate dot plot

```
DotPlot(test, features = feature)
```



## ## Use CellCycleScoring function to calculate cell cycle scores

```
marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
```



## ## scDC

```
library(scDC)
res_scDC_noClust <- scDC_noClustering(filtered_test$total.cluster, filtered_test$orig.ident, calCI = TRUE, calCI_method = c("BCa"),conf_level=0.95, ncores =6 )
barplotCI(res_scDC_noClust, c("control","CM"))

merged_data_subT <- merge(subcluster_T, test@meta.data, by='Barcode')
merged_data_subMono <- merge(subcluster_Mono, test@meta.data, by='Barcode')

#scDC of sub_T
subT_res_scDC_noClust <- scDC_noClustering(merged_data_subT$T, merged_data_subT$orig.ident, calCI = TRUE, calCI_method = c("BCa"),conf_level=0.95, ncores =6 )
barplotCI(subT_res_scDC_noClust, c("control","CM"))

#scDC of sub_Mono
subMono_res_scDC_noClust <- scDC_noClustering(merged_data_subMono$Mono, merged_data_subMono$orig.ident, calCI=TRUE,calCI_method=c("BCa"),conf_level=0.95, ncores =6 )
barplotCI(subMono_res_scDC_noClust, c("control","CM"))

```



## ## AddModuleScore 

```
group <- readxl::read_xlsx("/gene_list.xlsx")
gene <- as.list(group)
test <- AddModuleScore(object = test,features = gene,ctrl = 80,name = 'Group_Features')
colnames(test@meta.data)
test@meta.data$total.cluster <- factor(test@meta.data$total.cluster,
                                           levels = c(
"B cell","Monocyte","enrythroid cell","neutrophil","T cell","endothelial cell","epithelial cell","kupffer cell","Dendritic cell","Eosinophil cell","Hepatocytes","Plasmacytoid dendritic cell","unknown","dying cell"
))

VlnPlot(test,features= 'group score', group.by = "total.cluster", pt.size=0)
```




## ## Monocle2

#### #Read rds data

```
data <- readRDS(infile.rds)
```

#### #Build S4 object

```
HSMM <- newCellDataSet(count, phenoData = pd, featureData = fd)
```

#### #Estimate size factors and dispersions

```
HSMM <- estimateSizeFactors(HSMM) 
HSMM <- estimateDispersions(HSMM)
```

#### #Filter low-quality genes

```
HSMM <- detectGenes(HSMM, min_expr = mins)
```

#### #Select genes

```
diff_test_res <- differentialGeneTest(HSMM[expressed_genes, ], fullModelFormulaStr = "~ Clusters") 
HSMM <- setOrderingFilter(HSMM, ordering_gene)
```

#### #Dimensionality reduction and ordering

```
HSMM <- reduceDimension(HSMM, max_components = max_components, method = method) HSMM <- orderCells(HSMM)
```

#### #Visualization

```
plot_cell_trajectory(HSMM, color_by = "Cluster") 
plot_cell_trajectory(HSMM, color_by = "Sample") 
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
```



## ## RNA Velocity Analysis (Velocyto(0.17.17)+ScVelo(0.2.3))

#### #Create loom file using Velocyto

```
velocyto run10x -m $rmsk_gtf $cellranger_outDir $cellranger_gtf
```

#### #Import packages

```
import anndata 
import scvelo as scv
```

#### #Read file

```
adata = anndata.read_loom(loom)
```

#### #Prepare data

```
scv.settings.verbosity = 3   #show errors(0), warnings(1), info(2), hints(3) scv.settings.presenter_view = True  #set max width size for presenter view scv.settings.set_figure_params('scvelo')  #for beautified visualization 

adata = scv.datasets.pancreas() 

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) 

scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
```

#### #Build dynamical model

```
scv.tl.recover_dynamics(adata) 
scv.tl.velocity(adata, mode='dynamical') 
scv.tl.velocity_graph(adata)
```

#### #Visualization

```
scv.pl.velocity_embedding_stream(adata, basis='umap') 

top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300] 

scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='clusters', n_convolve=100) 

top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index 

scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, frameon=False)
```





