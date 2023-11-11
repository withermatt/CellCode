
main_dir <- "/Users/brittanymarcelo/Documents/GitHub/CellCode/"
setwd(main_dir)

#load packages and functions
source("packages_functions.R")
plot_cluster_composition_barchart <- function(df, Freq_minimum = 0.1, var1, var2) {
  
  ggData <- as.data.frame(prop.table(table(df[,var1], df[,var2]), margin = 2))
  
  top <- subset(ggData, Freq >= Freq_minimum)
  top$Var2 <- factor(top$Var2)
  top$Var1 <- factor(top$Var1)
  
  
  p <- ggplot(data = top, aes(x = Var1, y = Freq, fill = Var2)) +
    geom_bar(position="stack", stat = "identity", show.legend = T) +
    labs(x = NULL, y = "Fraction") +
    scale_y_continuous(breaks = seq(0,0.8,0.2)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  return(p)
}

plot_cluster_composition_heatmap <- function(df, column_var, row_var) {
  
  #This normalizes the frequencies for the row variable (i.e each row sums to 1)
  ggData <- t(as.matrix(prop.table(table(df[,column_var], df[,row_var]), margin = 2)))
  
  hm <- ComplexHeatmap::Heatmap(ggData,
                                column_title = column_var,
                                column_title_side = "bottom",
                                row_title = row_var,
                                row_title_side = "right",
                                column_names_gp = gpar(fontsize = 4),
                                name = "Freq")
  
  return(hm)
  
}
GO_enrichment_analysis <- function(gene_list, list_names, target_id_source = NULL, sim_cutoff = 0.7) {
  
  if (is.null(target_id_source) == FALSE) {
    ind <- match(gene_list, target_id_source$gene_short_name)
    gene_list_ENS <- target_id_source$id[ind]
  } else {
    gene_list_ENS <- gene_list
  }
  
  ego <- enrichGO(gene = gene_list_ENS,
                  OrgDb = org.Mm.eg.db,
                  ont = "BP",
                  keyType = "ENSEMBL",
                  readable = TRUE)
  
  ego2 <- simplify(ego, cutoff = sim_cutoff)
  ego2 <- as.data.frame(ego2)
  if (nrow(ego2) > 0) {
    ego2$GeneRatio <- DOSE::parse_ratio(ego2$GeneRatio)
    ego2 <- ego2[!is.na(ego2$Description), ]
    ego2$cluster <- list_names
  } else {
    ego2 <- NULL
  }
  return(ego2)
}
filter_GOdata <- function(df, n = 20, GeneRatio_thresh = 0) {
  
  if (!is.null(df)) {
    df <- subset(df, GeneRatio >= GeneRatio_thresh)
    df <- df[order(-df$GeneRatio),]
    if (length(df$GeneRatio) >= n) {
      df <- df[1:n, ]
    }
  }
  
  return(df)
  
}

cds_processed <- readRDS("Robjects/cds_processed.rds")

#Cluster UMAP
cds_processed <- cluster_cells(cds_processed, reduction_method = "UMAP", k = 700, cluster_method = "louvain")
colData(cds_processed)$cluster700 <- clusters(cds_processed)
plot_cells(cds_processed, alpha = 0.5, group_label_size = 8, color_cells_by = "cluster1000")
#Plot cluster and condition proportions
pData1 <- as.data.frame(colData(cds_processed))
plot_cluster_composition_heatmap(df = pData1, var1 = "cluster500", var2 = "unique")


umap_coords <- as.data.frame(reducedDims(cds_processed)$UMAP)
colData(cds_processed)$umap_1 <- umap_coords$V1
colData(cds_processed)$umap_2 <- umap_coords$V2
saveRDS(cds_processed, "Robjects/cds_processed_clusters.rds")


#Repeat the clustering for single base condition
cds_by_base <- readRDS("Robjects/cds_by_base.rds")
cds3 <- cds_by_base[[3]]
cds3 <- cluster_cells(cds3, reduction_method = "UMAP", k = 600, cluster_method = "louvain")
colData(cds3)$cluster <- clusters(cds3)
plot_cells(cds3, alpha = 0.5, group_label_size = 8)
pData3 <- as.data.frame(colData(cds3))
plot_cluster_composition_heatmap(df = pData3, var2 = "cytokine", var1 = "cluster")









markers <- top_markers(cds_processed, group_cells_by="cluster", marker_sig_test = T)

#Gene ontology
ENS_source <- as.data.frame(rowData(cds_processed))
marker_genes <- lapply(split(markers, markers$cell_group), function(x){ x$gene_id })

#Perform Gene ontology enrichment analysis
GOdata <- pbmapply(GO_enrichment_analysis, gene_list = marker_genes, list_names = 1:length(marker_genes), sim_cutoff = 0.5, SIMPLIFY = F)
#Filter results
GOdata <- lapply(GOdata, filter_GOdata, n = 2)

#saveRDS(GOdata, "Robjects/GO_gm_results.rds")
#GOdata <- readRDS("Robjects/GO_gm_results.rds")

GOdata2 <- bind_rows(GOdata)
GOdata2$cluster <- as.factor(GOdata2$cluster)
GOdata2$Description <- as.factor(GOdata2$Description)
GOdata2$qvalue2 <- log10(GOdata2$qvalue)*-1

#GOdata2 <- subset(GOdata2, Count >=5 | GeneRatio >= 0.4)
#GO_subset <- subset(GOdata2, cluster %in% c(2,7,26,42,35,23,15,3,20,41,39,14,12, 11, 8, 19))

ggplot(data = GOdata2, aes(x = GeneRatio, y = Description, size = Count, color = qvalue2)) +
  geom_point() +
  scale_color_continuous(low="black", high="red",
                         guide=guide_colorbar(reverse=F)) +
  facet_grid(cluster~., scale="free") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())







