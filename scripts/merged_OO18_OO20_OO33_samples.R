#clean up your environment
rm(list=ls())

install.packages('devtools')
devtools::install_github('immunogenomics/presto')

# Load required libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)
library(SingleR)
library(celldex)

# Load Seurat objects
seurat_object1 <- readRDS(file = "C:/Users/Administrator/Desktop/seurat analysis/OO18_sample/seurat_object_OO18_sample.rds")
seurat_object2 <- readRDS(file = "C:/Users/Administrator/Desktop/seurat analysis/OO20_sample/seurat_object_OO20_sample.rds")
seurat_object3 <- readRDS(file = "C:/Users/Administrator/Desktop/seurat analysis/OO33_sample/seurat_object_OO33_sample.rds")

# Check metadata of each Seurat object
head(seurat_object1@meta.data)
head(seurat_object2@meta.data)
head(seurat_object3@meta.data)

# Merge Seurat objects
seurat_object <- merge(seurat_object1, y = list(seurat_object2, seurat_object3), add.cell.ids = c("OO18_sample", "OO20_sample", "OO33_sample"), project = "merged")

# Check metadata of merged_seurat_object
head(seurat_object)

# Define the file path where you want to save the Seurat object
#save_path <- "C:/Users/Administrator/Desktop/seurat analysis/merged_OO18_OO20_OO33_samples/seurat_object_merged_OO18_OO20_OO33_samples.rds"  

# Save the Seurat object
saveRDS(seurat_object, file = "C:/Users/Administrator/Desktop/seurat analysis/merged_OO18_OO20_OO33_samples/seurat_object_merged_OO18_OO20_OO33_samples.rds")

seurat_object

head(seurat_object@meta.data)


# Load the saved Seurat object
seurat_object <- readRDS(file = "C:/Users/Administrator/Desktop/seurat analysis/merged_OO18_OO20_OO33_samples/seurat_object_merged_OO18_OO20_OO33_samples.rds")

head(seurat_object@meta.data)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(seurat_object@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# Filter cells that have unique feature counts over 2,500 or less than 200
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

# Filter cells that have >5% mitochondrial counts
seurat_object <- subset(seurat_object, subset = percent.mt < 5)

# Apply LogNormalize which includes a log transformation
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)


# Scaling the data
#all.genes <- rownames(seurat_object)
#seurat_object <- ScaleData(seurat_object, features = all.genes)
#seurat_object <- ScaleData(seurat_object, vars.to.regress = "percent.mt")

#scale all genes and regress out mitochondrial content simultaneously
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes, vars.to.regress = "percent.mt")

# Save the Seurat object
saveRDS(seurat_object, file = "C:/Users/Administrator/Desktop/seurat analysis/merged_OO18_OO20_OO33_samples/seurat_object_merged_OO18_OO20_OO33_samples_02.rds")

# Load the saved Seurat object
seurat_object <- readRDS(file = "C:/Users/Administrator/Desktop/seurat analysis/merged_OO18_OO20_OO33_samples/seurat_object_merged_OO18_OO20_OO33_samples_02.rds")

#Identification of highly variable features (feature selection)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

png("C:/Users/Administrator/Desktop/seurat analysis/merged_OO18_OO20_OO33_samples/plot_merged_variable_feature_2000.png", width = 1000, height = 800)
plot1 + plot2
dev.off()

# Perform linear dimensional reduction PCA
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
# Examine and visualize PCA results a few different ways
print(seurat_object[["pca"]], dims = 1:5, nfeatures = 5)

# visualize PCA
VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca")
DimPlot(seurat_object, reduction = "pca") + NoLegend()
DimPlot(seurat_object, reduction = "pca") 

DimHeatmap(seurat_object, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(seurat_object, dims = 1:15, cells = 500, balanced = TRUE)


png("C:/Users/Administrator/Desktop/seurat analysis/merged_OO18_OO20_OO33_samples/heatmap_merged_pca_1_to_6.png", width = 1600, height = 1600, res=150)
DimHeatmap(seurat_object, dims = 1:6, cells = 500, balanced = TRUE)

dev.off() 

# Determine the ‘dimensionality’ of the dataset
ElbowPlot(seurat_object)

# Cluster the cells
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat_object), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)
seurat_object <- RunUMAP(seurat_object, dims = 1:10)
seurat_object <- RunTSNE(seurat_object, dims = 1:10)

# visualize umap, tsne clusters
DimPlot(seurat_object, reduction = "umap")
DimPlot(seurat_object, reduction = "umap", label = TRUE) + NoLegend()
DimPlot(seurat_object, reduction = "tsne")
DimPlot(seurat_object, reduction = "tsne", label = TRUE) + NoLegend()

# Save the Seurat object
saveRDS(seurat_object, file = "C:/Users/Administrator/Desktop/seurat analysis/merged_OO18_OO20_OO33_samples/seurat_object_merged_OO18_OO20_OO33_samples_03.rds")

# Load the saved Seurat object
seurat_object <- readRDS(file = "C:/Users/Administrator/Desktop/seurat analysis/merged_OO18_OO20_OO33_samples/seurat_object_merged_OO18_OO20_OO33_samples_03.rds")

# Join layers for differential expression analysis
seurat_object <- JoinLayers(
  seurat_object, 
  layer.names = c("counts.OO18_sample", "counts.OO20_sample", "counts.OO33_sample"), 
  new.layer.name = "counts"
)

# Finding differentially expressed features (cluster biomarkers)
# find markers for every cluster compared to all remaining cells, report only the positive one
seurat_object.markers <- FindAllMarkers(seurat_object, only.pos = TRUE)
seurat_object.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Assuming your filtered data is still in seurat_object.markers and you want to print it directly:
seurat_object.markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  print(n = 20)  # This prints the first 20 rows of the filtered data.
head(seurat_object.markers)


# find all markers of cluster 0
cluster0.markers <- FindMarkers(seurat_object, ident.1 = 0)
cluster0.markers
head(cluster0.markers, n = 10)

# find all markers distinguishing cluster 0 from clusters 1 and 2
cluster0.markers <- FindMarkers(seurat_object, ident.1 = 0, ident.2 = c(1, 2))
head(cluster0.markers, n = 5)

# ROC test
cluster0.markers <- FindMarkers(seurat_object, ident.1 = 7, logfc.threshold = 2, test.use = "roc", only.pos = TRUE)
head(cluster0.markers, n = 10)
cluster0.markers

# visualizing marker expression
VlnPlot(seurat_object, features = c("CD14", "NOG", "CCR6", "LEF1", "GZMA", "ARHGAP15","NKG7", "LTB", "MS4A1","TMEM40"))

# plot raw counts as well
VlnPlot(seurat_object, features = c("CD14", "NOG", "CCR6", "LEF1", "GZMA", "ARHGAP15","NKG7", "LTB", "MS4A1","TMEM40"), layer = "counts", log = TRUE)

## visualizing feature expression on a umap 
FeaturePlot(seurat_object, features = c("CD14", "NOG", "CCR6", "LEF1", "GZMA", "ARHGAP15","NKG7", "LTB", "MS4A1","TMEM40"))

# expression heatmap for given cells and features
seurat_object.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
#DoHeatmap(seurat_object, features = top10$gene) + NoLegend()
DoHeatmap(seurat_object, features = top10$gene) + NoLegend() + theme(axis.text.y = element_text(size = 4))  # Adjust the size as needed
DoHeatmap(seurat_object, features = top10$gene) + theme(axis.text.y = element_text(size = 4))  # Adjust the size as needed


# annotation cluster
# Number of top markers to select per cluster
n_top_markers <- 50

# Arrange by p-value and avg_log2FC, then select the top markers for each cluster
top_markers_per_cluster <- seurat_object.markers %>%
  filter(avg_log2FC > 1) %>%
  arrange(cluster, p_val, desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  slice_head(n = n_top_markers) %>%
  ungroup()

# View the top markers for each cluster
print(top_markers_per_cluster)

# If you also want to view the top 20 rows regardless of cluster, just to get a broader overview
print(top_markers_per_cluster, n = 200)

# Identify the top 50 markers per cluster, just gene names
top10_genes_per_cluster <- top_markers_per_cluster %>%
  filter(avg_log2FC > 1) %>%
  arrange(cluster, p_val, desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  slice_head(n = 50) %>%
  summarise(genes = list(gene)) %>%
  ungroup()

# Using a loop to print each cluster's top 50 genes with correct cluster numbering starting from 0
for(i in 1:nrow(top10_genes_per_cluster)) {
  cat("Cluster", as.numeric(top10_genes_per_cluster$cluster[i]) - 1, ":")
  cat(paste(top10_genes_per_cluster$genes[[i]], collapse = ", "), "\n\n")
}


# visualizing marker expression
VlnPlot(seurat_object, features = c("FCAR", "RPS5"))

# plot raw counts as well
VlnPlot(seurat_object, features = c("FCAR", "RPS5"), layer = "counts", log = TRUE)

# visualizes feature expression on a umap to check top maker genes
FeaturePlot(seurat_object, features = c("FCAR", "RPS5", "IL7R", "BACH2", "CCL5", "AL136456.1", "LINGO2 ", "LTB","CCR10","TMEM40"))


## SingleR annotation
# Get normalized data for SingleR; ensure it has gene names as rownames
data <- GetAssayData(seurat_object, assay = "RNA", slot = "data")
if (is.null(rownames(data))) {
  rownames(data) <- rownames(seurat_object[["RNA"]])
}

# Load reference dataset
ref <- celldex::HumanPrimaryCellAtlasData()

# Perform SingleR cell annotation
annotations <- SingleR(test = as.SingleCellExperiment(seurat_object),ref = ref, labels = ref$label.main)

# Check the results
head(annotations)
annotations

# Adding annotations into seurat_object
# Convert the cell type predictions into a factor for better handling in Seurat
annotations$predicted.labels <- factor(annotations$labels)

# Add the annotations to the Seurat object's metadata
seurat_object$SingleR.labels <- annotations$predicted.labels[match(colnames(seurat_object), rownames(annotations))]

# check the added annotations
head(seurat_object$SingleR.labels)

# visualize the distribution of predicted cell types
DimPlot(seurat_object, group.by = "SingleR.labels")

# Visualize the distribution of predicted cell types with labels 
seurat_object <- SetIdent(seurat_object, value = "SingleR.labels")
DimPlot(seurat_object, label = TRUE, repel = TRUE, label.size = 4) + NoLegend()
DimPlot(seurat_object, label = TRUE, repel = TRUE, label.size = 4)

## Visualize the distribution of predicted cell types with cluster numbers
# Create a new identifier that combines Seurat cluster and SingleR labels
seurat_object$Combined.labels <- paste("Cluster", seurat_object$seurat_clusters, seurat_object$SingleR.labels, sep=": ")

# Set this new label as the active identity for plotting
seurat_object <- SetIdent(seurat_object, value = "Combined.labels")

# Visualize the distribution with the combined labels
plot <- DimPlot(seurat_object, label = TRUE, repel = TRUE, label.size = 3) + NoLegend()

# Print the plot
print(plot)

plot <- DimPlot(seurat_object, label = TRUE, repel = TRUE, label.size = 3) 

# Print the plot
print(plot)

# Create a new identifier that combines Seurat cluster numbers and SingleR labels directly
seurat_object$Combined.labels <- paste(seurat_object$seurat_clusters, seurat_object$SingleR.labels, sep=": ")

# Set this new label as the active identity for plotting
seurat_object <- SetIdent(seurat_object, value = "Combined.labels")

# Visualize the distribution with the combined labels
plot <- DimPlot(seurat_object, label = TRUE, repel = TRUE, label.size = 4) + NoLegend()

# Print the plot
print(plot)

plot <- DimPlot(seurat_object, label = TRUE, repel = TRUE, label.size = 4)

# Print the plot
print(plot)

# integrated_samples_SingleR_Cell_Type_Annotation_Across_Clusters
head(seurat_object@meta.data)
meta_data <- seurat_object@meta.data

# Ensure 'meta_data' is a dataframe and the columns 
print(colnames(meta_data))

# Create a data frame for plotting
data_to_plot <- data.frame(
  Cluster = meta_data$seurat_clusters,
  CellType = meta_data$SingleR.labels
)

# calculate proportions
plot_data <- data_to_plot %>%
  group_by(Cluster, CellType) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Cluster) %>%
  mutate(Proportion = Count / sum(Count))

# Check the resulting data frame
head(plot_data)

# Plot
ggplot(plot_data, aes(x = as.factor(Cluster), y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Percentage of Cell Types", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## manual annotation
# Set 'seurat_clusters' as the active identity
seurat_object <- SetIdent(seurat_object, value = "seurat_clusters")

# Define new cluster IDs and apply them
new.cluster.ids <- c("CD14+ Mono", "Naive CD4 T", "CD4 T", "Naive CD4 T", "NK_1 subset", "CD4 T", "NK_1 subset", "T",
                     "B", "Platelet")

names(new.cluster.ids) <- levels(seurat_object)
seurat_object <- RenameIdents(seurat_object, new.cluster.ids)

# Convert annotations to a factor and add them to Seurat object
annotations$predicted.labels <- factor(annotations$labels)
seurat_object[["new_annotations"]] <- annotations$predicted.labels[match(colnames(seurat_object), rownames(annotations))]

# Check the newly added annotations
head(seurat_object[["new_annotations"]])

# Check the updated metadata
head(seurat_object@meta.data, 5)

# Visualize the distribution of predicted cell types with labels
DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
plot <- DimPlot(seurat_object, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 18), legend.text = element_text(size = 15)) + guides(colour = guide_legend(override.aes = list(size = 5))) +
  xlim(c(-10, 15)) + ylim(c(-25, 15))
ggsave(filename = "C:/Users/Administrator/Desktop/seurat analysis/merged_OO18_OO20_OO33_samples/merged_manual_annotation_legend_label.jpg", height = 7, width = 12, plot = plot, quality = 50)


# DotPlot for genes expression in each cluster
genes <- c("CD14", "LYZ", "FCN1", "CCR7", "TCF7", "NOG", "IL7R", "CCR6", "LEF1", "GZMH", "GZMA", "GZMB", "CST7", "INPP4B", "ARHGAP15", "GNLY",
              "SPON2", "NKG7", "LTB", "TRAC", "CD27", "MS4A1", "TMEM40")

# Add spaces between each character in gene names
spaced_genes <- sapply(genes, function(gene) paste(unlist(strsplit(gene, "")), collapse = " "))

DotPlot(object = seurat_object, features = genes) +
theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 7))

# Create the DotPlot using these custom labels
DotPlot(seurat_object, features = genes) +
  scale_x_discrete(labels = spaced_genes) +  # Apply spaced gene names to x-axis
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 7))

# Save the Seurat object
saveRDS(seurat_object, file = "C:/Users/Administrator/Desktop/seurat analysis/merged_OO18_OO20_OO33_samples/seurat_object_merged_OO18_OO20_OO33_samples_04.rds")

# Load the saved Seurat object
seurat_object <- readRDS(file = "C:/Users/Administrator/Desktop/seurat analysis/merged_OO18_OO20_OO33_samples/seurat_object_merged_OO18_OO20_OO33_samples_04.rds")






