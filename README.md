## **Comparative scRNA-Seq Analysis Project Using Seurat**

### **Overview**

Single-cell RNA sequencing (scRNA-seq) has revolutionized our understanding of the complexities and dynamics of cellular processes at an unprecedented resolution. This project involves scRNA-seq data analysis of individual, **merged**, and **integrated** samples using the Seurat package. The objective is to compare cell type distributions and gene expression profiles across different processing methods to evaluate the impact of **merging versus integrating** datasets. This approach helps understand the influence of batch effects and the effectiveness of integration techniques in scRNA-seq analyses.

#### **Objectives**

-   **Comparative Analysis:** To systematically compare the outcomes of merged versus integrated scRNA-seq data analysis strategies. This comparison is crucial in assessing how different approaches influence the interpretation of cellular heterogeneity and the detection of distinct biological states.

-   **Batch Effect Evaluation:** To evaluate the effectiveness of data integration methods offered by Seurat in mitigating batch effects, which are known to skew analysis results and obscure genuine biological variations.

-   **Optimization of scRNA-Seq Workflows:** To refine analytical workflows in scRNA-seq studies, ensuring that the data processing techniques enhance the biological interpretability of the data without introducing artifacts.

#### **Methodological Approach**

-   **Data Processing:** Utilize Seurat, a widely acclaimed R package designed for QC, analysis, and exploration of scRNA-seq data. The project involves:

    -   **Data Cleaning and Normalization:** Standardizing data to remove technical noise.

    -   **Data Merging Process:** Combining datasets from individual samples into a single dataset. This process involves aligning data from different samples or conditions without the explicit modeling of batch effects, **which can be beneficial for increasing statistical power** but **may retain batch-related variations**.

    -   **Data Integration:** **Harmonizing datasets** from multiple batches or experiments using advanced algorithms to **minimize batch effects** and improve comparability. This is crucial for analyses where batch effects can obscure biological insights.

    -   **Clustering and Differential Expression Analysis:** Identifying and characterizing cell clusters based on their gene expression profiles to detect differentially expressed genes across conditions.

-   **Statistical Analysis:** Employ advanced statistical techniques to discern meaningful patterns and relationships in the data, facilitating a robust understanding of the underlying biological processes.

### **Datasets**

#### **Individual Samples:**

-   **OO18_sample:** Peripheral blood mononuclear cells (PBMC) from dataset OO18.

-   **OO20_sample:** Peripheral blood mononuclear cells (PBMC) from dataset OO20.

-   **OO33_sample:** Peripheral blood mononuclear cells (PBMC) from dataset OO33.

Each of these samples represents a unique collection of PBMCs, possibly under different conditions or from different individuals, providing a diverse basis for analyzing cellular responses and interactions.

#### **Merged Sample:**

-   **Merged PBMC Sample:** Combines OO18, OO20, and OO33 using Seurat's merge function. This merged dataset aggregates the individual datasets without adjusting for potential batch effects, offering a broader view of the combined cellular landscape across all samples.

#### **Integrated Sample:**

-   **Integrated PBMC Sample:** Integrates OO18, OO20, and OO33 using Seurat's data integration and batch correction methods. This approach is designed to minimize batch effects while preserving the genuine biological signals across the datasets, allowing for a more accurate comparison and analysis of cellular dynamics and gene expression patterns within PBMCs.

#### **Importance of PBMCs in Research:**

Peripheral blood mononuclear cells (PBMCs) are a critical component in immunological research due to their role in the immune system and their accessibility as a sample type. PBMCs include lymphocytes (T cells, B cells, and NK cells), monocytes, and dendritic cells. Studies involving PBMCs are fundamental for understanding immune responses under various health and disease conditions, making them invaluable for research in immunology, vaccine development, and therapeutic testing.

The use of PBMCs in this project allows for the exploration of cell-specific responses and interactions, providing insights that are vital for advancements in medical research and clinical applications. The comparative analysis of merged versus integrated datasets further enhances our understanding of how data processing techniques can influence the interpretation of such complex biological data.

### **Analysis Workflow**

#### **1. Initial Setup and Data Loading:**

-   **Environment Setup:** Begin by clearing the R environment and installing necessary libraries such as Seurat, dplyr, ggplot2, and patchwork.

-   **Data Import:** Load individual Seurat objects for each dataset (OO18, OO20, OO33). This step includes reading data from source files, creating Seurat objects, and performing initial data inspections to ensure integrity and readiness for preprocessing.

#### **2. Quality Control:**

-   **Mitochondrial Content Assessment:** Calculate the percentage of mitochondrial genes to evaluate cell viability. Cells exhibiting unusually high mitochondrial content are typically indicative of stress or cell death and are filtered out.

-   **Cell Filtering:** Implement filters to exclude cells based on extreme values in the number of detected genes and mitochondrial content. This step ensures the retention of high-quality cells suitable for reliable analysis.

#### **3. Normalization and Data Scaling:**

-   **Normalization Techniques:** Employ the LogNormalize method for merged samples and SCTransform for integrated samples. These methods stabilize variance across the dataset, making it conducive for accurate downstream analysis.

-   **Identification of Variable Features:** Detect high-variance genes within each dataset, focusing on features likely to contribute significantly to data variability.

-   **Data Scaling:** Regress out confounding sources of variation, such as mitochondrial content, during data scaling. This process helps in reducing bias and improving the quality of clustering and other multivariate analyses.

#### **4. Dimensional Reduction and Clustering:**

-   **Principal Component Analysis (PCA):** Conduct PCA to reduce dataset dimensionality while capturing the major axes of variation.

-   **Clustering Preparation:** Analyze PCA results to select an optimal number of dimensions for clustering.

-   **UMAP/t-SNE Implementation:** Utilize UMAP or t-SNE based on the PCA outputs to visualize and explore data clustering.

-   **Cluster Identification:** Apply the Louvain algorithm to detect distinct community structures within the data, identifying potentially novel or known cell types.

#### **5. Differential Expression and Marker Identification:**

-   **Marker Gene Detection:** Use Seurat\'s **`FindAllMarkers`** function to determine genes that uniquely define each cluster, facilitating the identification of biological markers.

-   **Visualization Techniques:** Employ violin plots to display gene expression distribution and feature plots to visualize spatial expression patterns of identified markers across clusters.

#### **6. Annotation and Cell Type Identification:**

-   **Automated Cell Type Annotation:** Leverage **SingleR** with reference datasets such as the Human Primary Cell Atlas for automated cell type identification.

-   **Manual Annotation:** Enhance automated annotations by manually refining cell type labels based on **canonical markers** and expression profiles specific to identified clusters.

#### **7. Data Visualization and Output:**

-   **Clustering Visualization:** Generate UMAP and t-SNE plots to visually assess the effects of data integration versus merging, particularly focusing on the mitigation of batch effects.

-   **Expression Visualization:** Create violin and feature scatter plots to examine and compare gene expression across different clusters.

-   **Heatmap Generation:** Produce heatmaps to illustrate the expression patterns of top differentially expressed genes, providing insights into biological functions and interactions.

-   **Elbow and Dot Plots:** Use elbow plots to determine the number of principal components and dot plots to display the proportion and intensity of marker gene expression within each cluster.

-   **Differential Expression Analysis:** Prepare detailed visualizations and tables to showcase genes significantly upregulated in specific clusters, aiding in the interpretation of cellular functions and states.

#### **8. Saving and Documentation:**

-   **Data Saving:** Continuously save Seurat objects and analytical outputs at critical stages to ensure reproducibility and facilitate data recovery.

-   **Comprehensive Documentation:** Document all analysis steps, code, and results in a detailed manner to support findings and provide deep insights into data quality, analytical challenges, and solutions.

### **Conclusion**

#### **Key Findings**

-   **Batch Effect Mitigation:** The comparison between merged and integrated samples demonstrated that integration techniques, such as those implemented in Seurat, are effective at reducing batch effects. This is evident from the more homogeneous clustering of cell types in the integrated samples compared to the distinct, batch-influenced clusters in the merged samples.

-   **Cell Type Consistency:** Integrated samples showed a more consistent annotation of cell types across different datasets. This suggests that the integration process not only corrects for technical variations but also enhances the biological interpretability of the data.

-   **Gene Expression Profiling:** Differential expression analysis revealed distinct expression profiles between clusters, identifying specific marker genes that define cellular identities. This highlights the utility of scRNA-seq in uncovering nuanced biological differences at the single-cell level.

### **Figures by each step workflow**

#### Visualize QC metrics

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/merged_samples%20_feature-feature_scatter_plot_QC_metrics_violin_plot.jpg?raw=true){width="434"}

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/integrated_samples%20_feature-feature_scatter_plot_QC_metrics_violin_plot.jpg?raw=true){width="467"}

#### Identify the 10 most highly variable genes

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/merged_integrated_plot_merged_variable_feature_2000.jpg?raw=true)

#### Visualize PCA

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/heatmap_merged_pca_1_to_6.png?raw=true){width="572"}

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/heatmap_integration_pca_1_to_6.png?raw=true){width="569"}

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/VizDimLoadings_merged_integrated_pca_+ElbowPlot.jpg?raw=true)

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/VizDimLoadings_merged_integrated_pca.jpg?raw=true){width="519"}

#### SingleR annotation

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/OO18_sample_SingleR_annotation_legend_label_02.jpeg?raw=true){width="581"}

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/OO20_sample_SingleR_annotation_legend_label_02.jpeg?raw=true){width="586"}

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/OO33_sample_SingleR_annotation_legend_label_02.jpeg?raw=true){width="585"}

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/SingleR_merged_integrated_annotation_cluster_umap_legend_label.jpg?raw=true){width="585"}

#### SingleR annotation & Cell Type Annotation Across Clusters

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/merged_SingleR_annotation_legend_label_02_and_Cell_Type_Annotation_Across_Clusters.jpg?raw=true){width="592"}

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/integration_SingleR_annotation_legend_label_02_and_Cell_Type_Annotation_Across_Clusters.jpg?raw=true){width="587"}

#### Expression heatmap for given cells and top 10 features each cluster

**merged_samples**

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/heatmap_merged_cells_and_features_top_10_markers_for_each_cluster.jpeg?raw=true)

**integrated_samples**

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/integrated_samples_top10_genes_heatmap.jpeg?raw=true)

#### Visualizing marker expression merged_samples

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/merged_samples_violin_visualizing_marker_expression_each_cluster.jpeg?raw=true)

#### Visualizing marker row counts merged_samples

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/merged_sample_violin_visualizing_marker_expression_raw_counts_each_cluster.jpeg?raw=true)

#### Visualizing marker expression integrated_samples

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/integrated_samples_violin_visualizing_marker_expression_each_cluster.jpeg?raw=true)

#### Visualizing marker row counts integrated_samples

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/integrated_sample_violin_visualizing_marker_expression_raw_counts_each_cluster.jpeg?raw=true)

#### Visualizing feature expression on a umap merged_samples

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/FeaturePlot_merged_top_markers_per_cluster.jpeg?raw=true)

#### Visualizing feature expression on a umap integrated_samples

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/FeaturePlot_integrated_top_markers_per_cluster.jpeg?raw=true)

#### merged sample canonical markers annotation and DotPlot

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/merged_sample_canonical_markers_annotation_and_DotPlot.jpg?raw=true){width="583"}

#### integration sample canonical markers annotation & DotPlot

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/integration_sample_canonical_markers_annotation_and_DotPlot.jpg?raw=true){width="587"}

#### manual annotation by canonical markers

![](https://github.com/chingyaousf/Comparative-scRNA-Seq-Analysis-Project-Using-Seurat/blob/main/plots/merged_integrated_annotation_cluster_umap_legend_label.jpg?raw=true)

## Blog:

<https://ssidmarine.wordpress.com/2024/05/13/comparative-scrna-seq-analysis-project-using-seurat/>
