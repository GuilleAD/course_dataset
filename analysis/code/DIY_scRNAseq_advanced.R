# Introduction ----
# the goal of this script is to explore several advanced topics in single cell analysis, including working with multi-modal data (in this case, CITE-seq), subclustering, pseudobulking, and GO analysis.
# To explore these topics, we need yet another example dataset.  
# We'll use unpublished data from my lab that draws on how the mucosal immune system responds to different kinds of infections.
# You can download this dataset from the course lecture page or the data page.

# Lecture 1 major topics Multimodal data, Subclustering, Pseudobulk DEG, GO analysis 
# Load the Libraries----
library(Seurat) # does all the heavy lifting in this script, and needed to work with Seurat Objects
library(tidyverse)
library(gprofiler2) # used to convert human cell cycle genes to mouse orthologs
library(scCustomize) # alternative and often improved visuals for Seurat graphing functions
library(presto) # Seurat implements presto to enhance speed of "FindMarkers" type functions
library(plotly) # we'll use it here to make interactive UMAPs
library(ggrepel) # helps with our plot labels on a crowded UMAP space


# Load the Seurat Object ----
# In this experiment we are looking at the mesenteric lymph node from mice that are controls or that have been infected by one of several pathogens/microbes. For expediency, the cells in this project have already have preliminary coarse labels and fine/granular cell type labels. Additionally, the total cells have been dramatically reduced to work within the confines of a single class and resource limitations. 
MLN_seurat <- LoadSeuratRds("MLN_subsampled_DIY.Rds")
#SaveSeuratRds(MLN_seurat,"MLN_subsampled_DIY.Rds")

#Take a look at the Seurat object - what metadata is present to give us information about the samples?
MLN_seurat
unique(MLN_seurat$InfectionStatus)
unique(MLN_seurat$Tissue)

# Take a look at the object and notice that we now have 2 different assays "RNA" and "ADT". In this script, one of the things we will briefly review is how to handle data containing two assays or modalities and how to leverage that information in our clustering and dimensional reduction (UMAPS)

# Quality Control check of the data ----
VlnPlot(MLN_seurat, "nCount_RNA", group.by = "orig.ident", pt.size = 0) + NoLegend()
VlnPlot(MLN_seurat, "nFeature_RNA", group.by = "orig.ident", pt.size = 0)+ NoLegend()
VlnPlot(MLN_seurat, "mt.prop", group.by = "orig.ident", pt.size = 0) + NoLegend()
VlnPlot(MLN_seurat, "prot.size", group.by = "orig.ident", pt.size = 0) + NoLegend()#protein size is similar to nCount_RNA and is simply the log10 transformed read count from the antibody-based sequencing for each cell

# Optional - Advanced Visualizations of QC data
rna_counts <-c() # define empty obj
rna_counts_sd <- c()
Idents(MLN_seurat) <- "orig.ident" # change the working identities to the name given to each unique sample
#Extract the mean and standard deviation for transcript count for each individual sample
for (ident in  levels(MLN_seurat@active.ident)) {
  rna_counts <- c(rna_counts, mean(MLN_seurat@meta.data$nCount_RNA[MLN_seurat@meta.data$orig.ident==ident]))}
for (ident in  levels(MLN_seurat@active.ident)) {
  rna_counts_sd <- c(rna_counts_sd, sd(MLN_seurat@meta.data$nCount_RNA[MLN_seurat@meta.data$orig.ident==ident]))}
#Collect the data into a data frame
rna_counts <- as.data.frame(rna_counts)
rna_counts$orig.ident <- levels(MLN_seurat@active.ident)
rna_counts$SD <- rna_counts_sd

# Create a dataframe with unique 'orig.ident' and corresponding 'InfectionStatus'
metadata <- MLN_seurat@meta.data
df <- unique(metadata[, c("orig.ident", "InfectionStatus")])
rna_counts <- left_join(rna_counts, df, by = "orig.ident")

#Define a pallette with the approproate number of colors (we have 5 Infection types)
colors <- c("#980339", "#379FFB", "#D6A206", "#027102", "#6313B7" )

rna_counts %>% arrange(InfectionStatus) %>%
  mutate(orig.ident = factor(orig.ident, orig.ident)) %>%
  ggplot( aes( orig.ident, rna_counts)) + 
  geom_col(aes(fill = InfectionStatus)) + 
  scale_fill_manual(values=colors) + 
  ylab("RNA Counts +/- SD") +
  geom_errorbar(aes(ymin = rna_counts - SD, ymax = rna_counts + SD, alpha = 0.5)) +
  theme(panel.background = element_rect(fill = "white", 
            colour = NA), panel.border = element_rect(fill = NA, 
            colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(linewidth = rel(0.5)), 
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "grey85", 
        colour = "grey20"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# Perform integration on both modalities separately ----

#Perform Clustering and Integration on the RNA assay. Notice that right now, our RNA assay only has a "counts" layer. This means it has not been normalized or scaled yet.
# This entire section took 3 minutes on my laptop to get to the UMAP
DefaultAssay(MLN_seurat) <- "RNA"
MLN_seurat[['RNA']] <-  JoinLayers(MLN_seurat[['RNA']])
MLN_seurat[['RNA']] <-  split(x = MLN_seurat[['RNA']], f = MLN_seurat$InfectionStatus)
MLN_seurat <- NormalizeData(MLN_seurat)
MLN_seurat <- FindVariableFeatures(MLN_seurat, selection.method = "vst", nfeatures = 3000) 

#Visualize Variable Features
VariableFeaturePlot(MLN_seurat, cols = c("black", "red"),
                    pt.size = 0.1, log = NULL, selection.method = NULL, raster = FALSE) + scale_y_log10()
#Scale the data 
all.genes <- rownames(MLN_seurat[["RNA"]])
MLN_seurat <- ScaleData(MLN_seurat, features = all.genes)
MLN_seurat <- RunPCA(MLN_seurat, verbose = FALSE, reduction.name = "unintegrated_RNA_pca")

# Now that you have your PCA result above, we need to figure out how many of our PCs (i.e. how much of our variance) to use for UMAP and clustering.
# There is no single best way to pick the number of PCs to use for UMAP, but let's look at two common approaches. 
# Approach 1: Determine which PCs cumulatively account for 90% of the total variance
ElbowPlot(MLN_seurat, reduction = "unintegrated_RNA_pca", ndims = 50)
# Determine percent of variation associated with each PC
pct <- MLN_seurat[["unintegrated_RNA_pca"]]@stdev / sum(MLN_seurat[["unintegrated_RNA_pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
cumu
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # Wow!  That's a lot of PCs we need to consider!

# Approach 2 to determine the PC at which point the percentage of change in PC drops to <0.1 
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2 #gives a very low number - usually this is the one that is more accurate
#Visualize the total percent variation captured by the cumulative PC chosen
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

# Now that we have a good idea of how many PCs to use, we can run the clustering and UMAP.
MLN_seurat <- FindNeighbors(MLN_seurat, dims = 1:30, reduction = "unintegrated_RNA_pca", graph.name = "rna_unint")
MLN_seurat <- FindClusters(MLN_seurat, resolution = 2, verbose = FALSE, graph.name = "rna_unint", cluster.name = "unintegrated_RNA_clusters", group.singletons = T)
MLN_seurat <- RunUMAP(MLN_seurat, dims = 1:30, reduction = "unintegrated_RNA_pca", reduction.name = "unintegrated_RNA_umap")
DimPlot(MLN_seurat, reduction = "unintegrated_RNA_umap", group.by = c("unintegrated_RNA_clusters"), combine = T, label.size = 4, label = TRUE) + NoLegend()
DimPlot(MLN_seurat, reduction = "unintegrated_RNA_umap", group.by = c("InfectionStatus"), combine = T)

# Integrate the RNA across the different samples. Remember that the section below only works because we have split our data by infection and we have set the default assay to "RNA". 
# If these parameters are changed for any reason, they would need to be set back.
options(future.globals.maxSize = 40 * 1024^4) #The amount of memory available to you will entirely depend on the machine/laptop you have.

MLN_seurat <- IntegrateLayers(
  object = MLN_seurat, method = RPCAIntegration, 
  orig.reduction = "unintegrated_RNA_pca", new.reduction = "integrated.rpca",
  verbose = FALSE)#This took 1-2 minutes on my laptop

# remember, it is best practice to repeat FindNeighbors and FindClusters after integration
MLN_seurat <- FindNeighbors(MLN_seurat, reduction = "integrated.rpca", dims = 1:30, graph.name = "RNA_int")
MLN_seurat <- FindClusters(MLN_seurat, resolution = 2, cluster.name = "rpca_integration_clusters",graph.name = "RNA_int",  group.singletons = T)
MLN_seurat <- RunUMAP(MLN_seurat, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rna.rpca")
DimPlot(MLN_seurat, reduction = "umap.rna.rpca", group.by = c("rpca_integration_clusters"), combine = T, label.size = 4, label = TRUE) + NoLegend()
DimPlot(MLN_seurat, reduction = "umap.rna.rpca", group.by = c("InfectionStatus"), combine = T)

# Let's take a look at a gene of interest
FeaturePlot(MLN_seurat, "Cd4", reduction = "umap.rna.rpca") + ggtitle("Cd4 expression") 

# Cluster the Protein/Antibody based data 
# In addition to "counts", the Protein/ADT assay also contains a "data" slot. This is because normalization has already been performed for you using the dsb package. If you would like to read more about noise in CITE-seq assays and how to handle that noise, please see https://www.nature.com/articles/s41467-022-29356-8

# begin by changing the default assay to the protein data
DefaultAssay(MLN_seurat) <- "ADT"
# define proteins to use in clustering (remove the isotypes from the clustering list)
rownames(MLN_seurat)
Isotype_labels <- grep("isotype", rownames(MLN_seurat), ignore.case = TRUE, value = TRUE)
prots = rownames(MLN_seurat)[!rownames(MLN_seurat) %in% Isotype_labels] #remove isotypes
MLN_seurat[['ADT']] <-  JoinLayers(MLN_seurat[['ADT']])
MLN_seurat[['ADT']] <-  split(x = MLN_seurat[['ADT']], f = MLN_seurat$InfectionStatus)
MLN_seurat <- ScaleData(MLN_seurat, assay = "ADT")
MLN_seurat <- RunPCA(MLN_seurat, features =prots, reduction.name = "unintegrated_pca_adt", reduction.key = "pcaADT_", verbose = FALSE)
DimPlot(MLN_seurat, reduction = "unintegrated_pca_adt", group.by = "orig.ident")
ElbowPlot(MLN_seurat, reduction = "unintegrated_pca_adt", ndims = 40)

MLN_seurat <- FindNeighbors(object = MLN_seurat, dims = 1:20, features = prots, graph.name = "ADT_unint", 
                         k.param = 30, verbose = FALSE, reduction = "unintegrated_pca_adt" )

MLN_seurat <- FindClusters(MLN_seurat, resolution = 2,  graph.name = "ADT_unint", verbose = FALSE, cluster.name = "unintegrated_ADT_clusters", group.singletons = T)
MLN_seurat[['ADT']] <-  JoinLayers(MLN_seurat[['ADT']])
MLN_seurat <- RunUMAP(object = MLN_seurat, assay = "ADT", features = prots, 
                   n.neighbors = 30,  verbose = FALSE, reduction.name = "umap.ADT.unintegrated")

DimPlot(MLN_seurat, reduction = "umap.ADT.unintegrated", group.by = c("unintegrated_ADT_clusters", "InfectionStatus"),
        combine = FALSE, label.size = 2, label = TRUE)

# Integrate the ADT data 
MLN_seurat[["ADT"]] <- split(MLN_seurat[["ADT"]], f = MLN_seurat$InfectionStatus)
MLN_seurat <- IntegrateLayers(
  object = MLN_seurat, method = RPCAIntegration, features = prots, 
  orig.reduction = "unintegrated_pca_adt", new.reduction = "integrated.ADT.rpca", verbose = TRUE)
MLN_seurat <- FindNeighbors(MLN_seurat, reduction = "integrated.ADT.rpca", dims = 1:20,  graph.name = "ADT_int"  )
MLN_seurat <- FindClusters(MLN_seurat, resolution = 2, cluster.name = "adt_integration_clusters", graph.name = "ADT_int",  group.singletons = T)
MLN_seurat <- RunUMAP(MLN_seurat, reduction = "integrated.ADT.rpca", dims = 1:20, reduction.name = "umap.integrated.ADT")
DimPlot(MLN_seurat, reduction = "umap.integrated.ADT", group.by = c("adt_integration_clusters", "InfectionStatus"),combine = FALSE, label.size = 2)

# Just as we looked at a select gene above, we can now view individual protein expression across all the cells in our integrated datset
FeaturePlot(MLN_seurat, "CD4", reduction = "umap.integrated.ADT") + ggtitle("CD4 Protein") 
FeaturePlot(MLN_seurat, "CD8a", reduction = "umap.integrated.ADT") + ggtitle("CD8 Protein") 



# Combined (multimodal) analysis using Weighted Nearest Neighbor ----
# Note that this is the same approach you would use for analysis of multiome data (i.e. ATAC + RNA) or any other multi-modal data.
# You can learn about this method from this paper: https://doi.org/10.1016/j.cell.2021.04.048
# This section takes ~2-3 minutes 
DefaultAssay(MLN_seurat) <- "RNA"
MLN_seurat = FindMultiModalNeighbors(
  MLN_seurat, reduction.list = list("integrated.rpca", "integrated.ADT.rpca"), 
  dims.list = list(1:30, 1:20), 
  modality.weight.name = "RNA.weight", 
  knn.graph.name = "wknn",
  snn.graph.name = "wsnn",
  weighted.nn.name = "weighted.nn",
  verbose = FALSE)

MLN_seurat <- RunUMAP(MLN_seurat, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
MLN_seurat <- FindClusters(MLN_seurat, graph.name = "wsnn",  resolution = 2, verbose = FALSE, cluster.name = "WNN_clusters", n.start = 25, group.singletons = T)

DimPlot(MLN_seurat, reduction = "wnn.umap",   group.by = c( "WNN_clusters", "InfectionStatus"),
        combine = FALSE, label.size = 2, label = TRUE)

# Cell Cycle Scoring ----
# A list of cell cycle markers, from Kowalczyk et al, 2015, is loaded with Seurat (updated 2019).  
# We can segregate this list into markers of G2/M phase and markers of S phase

#Get the list of mouse cell cycle genes
m.s.genes = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
m.g2m.genes = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

DefaultAssay(MLN_seurat) <- "RNA"
MLN_seurat <- JoinLayers(MLN_seurat)
MLN_seurat <- CellCycleScoring(MLN_seurat, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
DimPlot(MLN_seurat, reduction = "wnn.umap", group.by = c("Phase"),combine = FALSE, label.size = 2)

# After further analysis, you might decide to regress out the cell cycle signature if it overwhelms cell type specific signatures. We will not do that here for the sake of time. 

# Additional single cell visualizations and quantification ----
# The scCustomize package has very similar functions to seurat but with different default settings 

FeaturePlot(MLN_seurat, "Cd4", reduction = "wnn.umap") + ggtitle("Cd4 expression") 
FeaturePlot_scCustom(MLN_seurat, "Cd4", reduction = "wnn.umap") + ggtitle("Cd4 expression") 
FeaturePlot_scCustom(MLN_seurat, "CD4", reduction = "wnn.umap") + ggtitle("CD4 Protein") 

DimPlot(MLN_seurat, reduction = "wnn.umap", group.by = c("WNN_clusters", "InfectionStatus"),combine = FALSE, label.size = 2)
DimPlot_scCustom(MLN_seurat, reduction = "wnn.umap", group.by = c("WNN_clusters", "InfectionStatus"),combine = FALSE, label.size = 2)


#Enhancing UMAPs with ggplot
# Get UMAP coordinates and cell types
umap_data <- as.data.frame(MLN_seurat[["wnn.umap"]]@cell.embeddings)
umap_data$cell_id <- rownames(umap_data)  
umap_data$cell_type <- MLN_seurat$CoarseCellType 
umap_data$FineLabel <- MLN_seurat$FineLabel 

# Create the ggplot
umap_plot <- ggplot(umap_data, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = cell_type, text = FineLabel)) +
  geom_point(size = 1.5) +
  theme_minimal() +
  labs(title = "UMAP of Cell Types") +
  theme(legend.position = "right")

print(umap_plot)

#plotly can make interactive versions of plots. Here we can toggle over the coarse annotation to see the fine annotation. 
interactive_umap <- ggplotly(umap_plot, tooltip = "text")
interactive_umap #Display the interactive plot

#Visualizing the UMAPS
umap_data <- MLN_seurat[["wnn.umap"]]@cell.embeddings
cell_type_mapping <- MLN_seurat$CoarseCellType

umap_data <- as.data.frame(umap_data)
umap_data$cluster <- factor(cell_type_mapping)

centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(wnnUMAP_1 = mean(wnnUMAP_1), wnnUMAP_2 = mean(wnnUMAP_2))


ggplot(umap_data, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.5) +  # Plot the points (adjust size & transparency)
  geom_text_repel(data = centroids, aes(label = cluster), size = 3, fontface = "bold", 
                  color = "black",
    box.padding = 0.1,  # Increase space around the label itself
    #point.padding = 10,  # Increase space between the label and data point
                  min.segment.length = 0.5) +
  scale_color_viridis_d() +  # Customize color scale
  theme_minimal() +  # Minimal theme for better presentation
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) 


# Cell type frequencies across infections
# Extract the relevant metadata from Seurat object
umap_data <- as.data.frame(MLN_seurat[["wnn.umap"]]@cell.embeddings)
umap_data$cluster <- factor(MLN_seurat$CoarseCellType)
umap_data$condition <- MLN_seurat$InfectionStatus  

# Group by 'condition' and 'cluster' and count the occurrences
cluster_condition_counts <- umap_data %>%
  group_by(condition, cluster) %>%
  summarise(count = n(), .groups = 'drop')

# Calculate the total number of cells in each condition
total_cells_per_condition <- umap_data %>%
  group_by(condition) %>%
  summarise(total_cells = n(), .groups = 'drop')

# Merge the counts with the total cells per condition
cluster_condition_percent <- cluster_condition_counts %>%
  left_join(total_cells_per_condition, by = "condition") %>%
  mutate(percentage = (count / total_cells) * 100)  # Calculate percentage

unique(MLN_seurat$InfectionStatus)
cluster_condition_percent$condition <- factor(cluster_condition_percent$condition, levels = c( "Naive" , "Candida", "Cryptosporidium", "H.poly",   "MNV",  "Yersinia" ))

# Create a stacked barplot of percent frequency of each cell type within each condition
plot <- ggplot(cluster_condition_percent, aes(x = cluster, y = percentage, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Condition", y = "Frequency") + 
  #scale_y_break(c(20, 60)) +
  scale_fill_manual(values = c("#F6A800", "#0072B2", "#D55E00", "#009E73", "#CC79A7")) +
  theme_minimal(base_size = 15) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.y = element_blank(),
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 12),  # Adjust legend text size
    legend.position = "bottom"  # Place the legend at the bottom
  ) + coord_flip()

plot


# Cell subsetting and reclustering-----
# Frequently with large data sets, the initial dimensional reductions and clustering is incapable of capturing the granularity provided by the data. 
# In this case, we can take a subset of interest in the data and recluster it. 
DimPlot(MLN_seurat, reduction = "wnn.umap", group.by = c("WNN_clusters", "Phase", "CoarseCellType"), combine = F)

#Subset clusters or a population of interest 
Idents(MLN_seurat) <- "CoarseCellType"
unique(Idents(MLN_seurat))
options(future.globals.maxSize = 40 * 1024^4) 
MLN_T <- subset(x = MLN_seurat, idents = c("T Cells and ILCs"))
# it is common practice to recluster & reintegrate data after subsetting
# we will now do this separately for the RNA and ADT assays.
# let's start with the RNA
# we begin with joinLayers because that ensures that we don't accidentally try to run split on data that is already split 
MLN_T[['RNA']] <-  JoinLayers(MLN_T[['RNA']])
MLN_T[['RNA']] <- split(x = MLN_T[['RNA']], f = MLN_T$InfectionStatus)
MLN_T <- FindVariableFeatures(MLN_T, selection.method = "vst", nfeatures = 3000) # Find new variable features across the selected populations
MLN_T <- RunPCA(MLN_T, verbose = FALSE, reduction.name = "T_RNA_pca", assay = "RNA" )
MLN_T <- IntegrateLayers( object = MLN_T, method = RPCAIntegration,  orig.reduction = "T_RNA_pca", new.reduction = "integrated.rpca.T", verbose = FALSE)
# now repeat this process for the ADT data
DefaultAssay(MLN_T) <- "ADT"
MLN_T[['ADT']] <-  JoinLayers(MLN_T[['ADT']])
MLN_T[['ADT']] <-  split(x = MLN_T[['ADT']], f = MLN_T$InfectionStatus)
Isotype_labels <- grep("isotype", rownames(MLN_T), ignore.case = TRUE, value = TRUE)
prots = rownames(MLN_T)[!rownames(MLN_T) %in% Isotype_labels] #remove isotypes
MLN_T <- RunPCA(MLN_T, features =prots, reduction.name = "T_ADT_pca", reduction.key = "pcaADT_", verbose = FALSE, assay = "ADT")
MLN_T <- IntegrateLayers( object = MLN_T, method = RPCAIntegration, features = prots,  orig.reduction = "T_ADT_pca", new.reduction = "integrated.ADT.rpca.T", verbose = TRUE)
# now we recluster the combined data
DefaultAssay(MLN_T) <- "RNA"
MLN_T = FindMultiModalNeighbors(
  MLN_T, reduction.list = list("integrated.rpca.T", "integrated.ADT.rpca.T"),  dims.list = list(1:25, 1:20),  modality.weight.name = "RNA.weight",  k.nn = 15, 
  knn.graph.name = "wknn_T",  snn.graph.name = "wsnn_T", weighted.nn.name = "weighted_T.nn", verbose = TRUE)

MLN_T <- RunUMAP(MLN_T, nn.name = "weighted_T.nn", reduction.name = "wnn.T.umap", reduction.key = "wnnUMAP", seed.use = 12, n.neighbors = 15L)
MLN_T <- FindClusters(MLN_T, graph.name = "wsnn_T",  resolution = 2, verbose = FALSE, cluster.name = "Tcell_clusters", n.start = 20, group.singletons = T)

#Visualize the results
DimPlot(MLN_T, reduction = "wnn.T.umap",   group.by = c( "Tcell_clusters", "FineLabel"),
        combine = FALSE, label.size = 2, label = TRUE)
FeaturePlot_scCustom(MLN_T, reduction = "wnn.T.umap", c("rb.prop", "nFeature_RNA")) # By reclustering the data, there are clearly a couple small clusters which are technically distinct and have very low ribosomal reads relative to the rest of the data. There's no right answer but in many cases, we may choose to remove these cells if we think they are technical artifacts
FeaturePlot_scCustom(MLN_T, reduction = "wnn.T.umap", c("Cxcr5", "Pdcd1", "Foxp3", "Rorc"))
FeaturePlot_scCustom(MLN_T, reduction = "wnn.T.umap", c("Eomes", "Tbx21", "Cd8a", "Cxcr6"))


# Pseudobulk Differential Expression ----
# In the previous lectures we used FindMarkers function to perform a baseline differential expression test across populations of interest. 
# However, these tests are known to greatly inflate statistical significance and when possible, pseudobulk analyses provide more accurate comparisons. 
# We will be doing this analysis with a minimum of 2 replicates per condition, but more replicates are preferred.

# To being, we pseudobulk the counts based on infection, mouse, and cell type
pseudo <- AggregateExpression(MLN_seurat, assays = "RNA", return.seurat = T, group.by = c("InfectionStatus", "Mouse", "FineLabel"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo))
pseudo$celltype.condition <- paste(pseudo$FineLabel, pseudo$InfectionStatus, sep = "_")

MLN_seurat$celltype.condition <- paste(MLN_seurat$FineLabel, MLN_seurat$InfectionStatus, sep = "_")
Idents(MLN_seurat) <-"celltype.condition"
mono_deg <- FindMarkers(MLN_seurat, ident.1 = "Monocyte Macrophage_Yersinia", 
                        ident.2 = "Monocyte Macrophage_Naive" )

Idents(pseudo) <- "celltype.condition"

bulk.mono.de <- FindMarkers(object = pseudo, 
                            ident.1 = "Monocyte Macrophage_Yersinia", 
                            ident.2 = "Monocyte Macrophage_Naive",
                            test.use = "DESeq2", 
                           min.cells.group = 2)
head(bulk.mono.de, n = 15)

sum(mono_deg$p_val_adj < 0.01)
sum(bulk.mono.de$p_val_adj < 0.01)




