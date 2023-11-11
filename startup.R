#Project:
  #CellCode in vitro screen

#Packages
library(monocle3)
#Tutorial for monocle3
#https://cole-trapnell-lab.github.io/monocle3/docs/introduction/

library(DescTools)
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
suppressPackageStartupMessages(library(ComplexHeatmap))
#Refernce manual for Heatmaps - very very helpful
#https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#titles-for-splitting

library(viridis)
library(circlize)
#library(randomForest)
library(caret)
library(pbmcapply)
library(pbapply)
library(patchwork)
#library(VennDiagram)
library(clusterProfiler)
library("org.Mm.eg.db", character.only = TRUE)
library(ClassDiscovery)
library(Matrix)
library(apcluster)
library(ggrepel)
library(glmnet)
library(dendextend)
library(clusterSim)
#library(factoextra)
library(lsa)
library(ggridges)


#User-defined functions

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summarized
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  #library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  if (!is.null(conf.interval)) {
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
  }
  
  return(datac)
}

compute_percentage <- function(x) {x/sum(x) * 100}

compute_fold_change <- function(x, standard) {log2(x/x[standard])}

RowRange <- function(x) {max(x) - min(x)}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#Need to update - may not be useful anymore
sc_pseudobulk <- function(cds, gene_groups, zscore = T, cell_group_column = "unique") {
  
  pData <- as.data.frame(colData(cds))
  rownames(pData(cds)) <- pData$barcode
  pData <- as.data.frame(colData(cds))
  
  cell_groups <- pData[,c("barcode", cell_group_column)]
  
  if (!is.null(gene_groups)) {
    agg0 <- aggregate_gene_expression(cds, cell_group_df = cell_groups, gene_group_df = gene_groups, norm_method = "log", scale_agg_values = zscore)
    
  } else {
    agg0 <- aggregate_gene_expression(cds, cell_group_df = cell_groups, norm_method = "log", scale_agg_values = zscore)
    
    #Swap target id for gene short name
    id2g <- as.data.frame(rowData(cds))
    old_row_names <- rownames(agg0)
    new_row_names <- match(old_row_names, id2g$id)
    rownames(agg0) <- id2g$gene_short_name[new_row_names]
    
    #agg0 <- agg0[gene_module_df$gene,]
    
  }
  
  return(as.matrix(agg0))
}

#Need to update
make_gene_module_hm <- function(mat, gene_module_df = NULL, aggregate_genes = T, title) {
  
  
  if (aggregate_genes) {
    
    hm <- Heatmap(mat,
                  show_row_dend = F,
                  show_column_dend = F,
                  show_column_names = T,
                  column_names_gp = gpar(fontsize = 5),
                  show_row_names = T,
                  row_names_side = "left",
                  cluster_columns = T,
                  cluster_rows = T,
                  border = T,
                  column_title = title,
                  row_title = "Gene Module",
                  name = "Expr")
    
  } else {
    
    split <- as.data.frame(gene_module_df$module)
    rownames(split) <- rownames(gene_module_df)
    colnames(split) <- "split"
    split$split <- factor(split$split, levels = c(1:max(gene_module_df$module)))
    
    hm <- Heatmap(mat,
                  show_row_dend = F,
                  show_column_dend = F,
                  show_column_names = T,
                  column_names_gp = gpar(fontsize = 5),
                  show_row_names = F,
                  row_names_side = "left",
                  cluster_columns = T,
                  cluster_rows = T,
                  #column_split = c(rep("Stim+IL-12", 27), rep("Stim", 29)),
                  #cluster_row_slices = F,
                  row_split = split,
                  row_gap = unit(3, "mm"),
                  border = T,
                  column_title = title,
                  #row_title = "Gene Module",
                  name = "Expr")
    
  }
  
  return(hm)
  
}

#Need to update
plot_aggregated_gene_modules <- function(mat, range_threshold = NULL, apcluster_pref = 0.3, return_grouping = F, normalize = F) {
  
  if(normalize) {
    agg_norm <- matrix(apply(mat, 1, compute_fold_change, standard = "Null_"), ncol = ncol(mat), nrow = nrow(mat), dimnames = dimnames(mat), byrow = T)
    gm_range <- apply(agg_norm, 1, RowRange)
    agg_norm2 <- agg_norm[,colnames(agg_norm) != "Null_"]
    if (!is.null(range_threshold)) {
      #Remove modules with minimal change
      agg_norm2 <- agg_norm2[unname(which(gm_range >= range_threshold)),]
    }
  } else {
    agg_norm2 <- mat
  }
  
  
  s <- negDistMat(t(agg_norm2))
  #Run clustering
  ap_cyt <- apcluster(s = s, q = apcluster_pref)
  groups <- ap_cyt@clusters
  ap_cyt2 <- aggExCluster(s = s, x = ap_cyt, includeSim = T)
  
  cyt_order <- names(unlist(groups[ap_cyt2@order]))
  
  #Create data frame to store cytokine groups
  cyt_groups <- bind_rows(mapply(FUN = function(x, i) {
    out_df <- data.frame(index = x, group = i, cytokine = names(x))
    return(out_df)
    
  }, x = groups, i = 1:length(groups), SIMPLIFY = F))
  
  #reorder matrix columns
  agg_norm3 <- agg_norm2[,cyt_order]
  
  #Add bar annotations to identify cytokine groups
  #Sort data frame to match matrix order
  cyt_groups_ordered <- cyt_groups[cyt_order,]
  cyt_groups_ha <- HeatmapAnnotation(Module = cyt_groups_ordered$group,
                                     col = list(Module = colorRamp2(1:length(groups), gg_color_hue(length(groups)))),
                                     simple_anno_size = unit(2, "mm"),
                                     show_annotation_name = FALSE, show_legend = F)
  
  col_split <- data.frame(split = factor(cyt_groups_ordered$group))
  
  hm <- Heatmap(agg_norm3,
                #col = colorRamp2(c(-1,-0.5,0,0.5,1), c("darkblue", "blue", "white", "red", "darkred")),
                show_row_dend = F,
                show_column_dend = F,
                show_column_names = T,
                column_names_gp = gpar(fontsize = 5),
                row_names_gp = gpar(fontsize = 8),
                show_row_names = T,
                row_names_side = "left",
                cluster_columns = F,
                cluster_rows = T,
                bottom_annotation = cyt_groups_ha,
                #right_annotation = rowAnnotation(FCrange = anno_barplot(gm_range), annotation_name_rot = 90, gap = unit(2, "mm")),
                column_split = col_split,
                column_title_side = "bottom",
                border = T,
                row_title = "Gene Module",
                name = "z-score")
  
  if (return_grouping) {
    return(list(hm, cyt_groups_ordered))
  } else {
    return(hm)
  }
  
}

#Convert monocle cds into Seurat object for use with Seurat functions
#I currently don't have a method to convert back. The easiest way is to add the new colData columns to colData of the original cds.
monocle_to_seurat <- function(monocle_cds) {
  requireNamespace("Seurat")
  data <- exprs(monocle_cds)
  export_cds <- Seurat::CreateSeuratObject(counts = data, 
                                           normalization.method = "LogNormalize",
                                           do.scale = TRUE, 
                                           do.center = TRUE,
                                           project = "exportCDS")
  export_cds@meta.data <- as.data.frame(colData(monocle_cds))
  
  return(export_cds)
}

## Plot heatmap of the composition of cells from each condition in each of the UMAP clusters.
##   df: a data frame.
##   column_var: defaults to cluster - must have this variable in the data frame to work (i.e run cluster_cells)
##   row_var: column indicating condition to group cells by. Rows get normalized to 1.
##   annotation_bar_groupings: column of data frame or vector of groupings for making a color bar annotation.
##   annotation_colors: color key for annotation bar. Must match length of unique group IDs from annotation_bar_groupings
##   row_label_annotation_thresh: maximum column (i.e cluster) enrichment needed for a row (i.e condition) to be annotated. 1 (default) hides all row names.
plot_cluster_composition_heatmap <- function(df, 
                                             column_var = "cluster", 
                                             row_var,
                                             cluster_columns = TRUE,
                                             column_reorder = TRUE,
                                             annotation_bar_groupings = NULL,
                                             annotation_colors = NULL,
                                             label_size = 6,
                                             row_label_annotation_thresh = 1) {
  
  #This normalizes the frequencies for the row variable (i.e each row sums to 1)
  hmData <- t(as.matrix(prop.table(table(df[,column_var], df[,row_var]), margin = 2)))
  hmData <- matrix(hmData, ncol = ncol(hmData), dimnames = dimnames(hmData))
  
  #Cluster columns
  
  # d <- hclust(dist(hmData), method = "average")
  # set.seed(12)
  # d <- reorder(as.dendrogram(d), hmData[,column_reorder_var])
  
  if (!is.logical(column_reorder)) {
    column_reorder <- hmData[column_reorder,] #which row to use for ordering weights
    
    #HACK to move another branch
    column_reorder[13] <- 1#column_reorder[2] * 1.5
    
  }
  
  #Determine which rows to annotate
  row_maximums <- rowMax(hmData)
  rows_to_label_pos <- which(row_maximums >= row_label_annotation_thresh)
  
  #annotate null conditions
  null_conditions <- c(unique(grep("Null", df$unique, value = TRUE)), unique(grep("CD3/28_IL-12_1ng/mL", df$unique, value = TRUE)))
  null_pos <- match(null_conditions, rownames(hmData))
  
  rows_to_label_pos <- unique(c(rows_to_label_pos, null_pos))
  rows_to_label <- rownames(hmData)[rows_to_label_pos]
  
  #Add annotation color bar
  if(!is.null(annotation_bar_groupings)) {
    
    if (length(annotation_bar_groupings) == 1) {
      classifier <- df[,annotation_bar_groupings][match(rownames(hmData), df[,row_var])]
    } else {
      classifier <- annotation_bar_groupings
    }
    
    #Set color key if none provided
    if (is.null(annotation_colors)) {
      annotation_colors = colorRamp2(1:length(unique(classifier)), sample(gg_color_hue(length(unique(classifier)))))
    }
    
    #Don't annotate specific rows for now
    # row_annot <- rowAnnotation(Group = factor(classifier),
    #                            Cond = anno_mark(at= rows_to_label_pos, labels = rows_to_label, which="row", side = "right", labels_gp = gpar(col= "black", fontsize = label_size)),
    #                             col = list(Group = annotation_colors),
    #                             simple_anno_size = unit(2, "mm"),
    #                             show_annotation_name = F, show_legend = T)
    
    row_annot <- rowAnnotation(Group = factor(classifier),
                               col = list(Group = annotation_colors),
                               simple_anno_size = unit(2, "mm"),
                               show_annotation_name = F, show_legend = T)
    
    
    hm <- ComplexHeatmap::Heatmap(hmData,
                                  column_title = column_var,
                                  column_title_side = "bottom",
                                  cluster_columns = cluster_columns,
                                  row_title = row_var,
                                  row_title_side = "right",
                                  show_row_names = TRUE,
                                  row_dend_reorder = TRUE,
                                  column_dend_reorder = column_reorder,
                                  column_names_gp = gpar(fontsize = 8),
                                  row_names_gp = gpar(fontsize = label_size),
                                  right_annotation = row_annot,
                                  name = "Freq")
    
  } else {
    
    row_annot <- rowAnnotation(Cond = anno_mark(at= rows_to_label_pos, labels = rows_to_label, which="row", side = "right", labels_gp = gpar(col= "black", fontsize = label_size)),
                               show_annotation_name = F, show_legend = T)
    
    hm <- ComplexHeatmap::Heatmap(hmData,
                                  column_title = column_var,
                                  column_title_side = "bottom",
                                  cluster_columns = cluster_columns,
                                  show_row_names = F,
                                  row_title = row_var,
                                  row_title_side = "right",
                                  row_dend_reorder = TRUE,
                                  column_dend_reorder = column_reorder,
                                  column_names_gp = gpar(fontsize = 8),
                                  #row_names_gp = gpar(fontsize = label_size),
                                  right_annotation = row_annot,
                                  name = "Freq")
  }
  
  return(hm)
  
}


## Fisher exact test for cluster enrichment by cells from a single condition
##   df: a data frame.
##   group_cells_by: column variable name indicating how to group cells.
##   split_by: Column variable to split the cell groups by. Each subgroup is tested against its own control. If NULL (default), every condition is tested against the same control condition (named by control_identifier).
##   control_identifier: If split_by is used, a string matching to a unique condition(s) to test against. If split_by is not used, a string matching exactly to one and only one group_cells_by names.
##   enrichment_probs: if TRUE (default), test for enrichment. If FALSE, test probabilities less than the null distribution to determine cluster depletion instead of enrichment.
test_cluster_enrichment <- function(df, 
                                    group_cells_by = "unique", 
                                    split_by = NULL, 
                                    control_identifier,
                                    greater_than = TRUE) {
  
  if (!is.null(split_by)) {
    dfs <- split(df, df[,split_by])
    
    results2 <- bind_rows(lapply(dfs, FUN = function(df_x) {
      
      ctrl <- unique(grep(control_identifier, df_x[,group_cells_by], value = TRUE))
      
      #Create counts table of clustering assignments
      counts_true <- as.data.frame(table(df_x[,"cluster"], df_x[,group_cells_by]))
      counts_split <- split(counts_true, counts_true$Var2)
      
      #Iterate through each condition to get a p-value from Fisher's exact test for significant enrichment 
      #of cells from that condition in each cluster relative to the control condition
      results <- pblapply(counts_split, FUN = function(x_test, x_control) {
        
        #Initialize dataframe to store probabilities for each cluster
        probs <- data.frame()
        
        x_total <- bind_rows(x_test,x_control)
        
        #Now iterate over each cluster and get probabilities and p-values
        for (clust in levels(x_test$Var1)) {
          
          # Initialize variables
          m <- sum(x_total$Freq[which(x_total$Var1 == clust)])    # cells IN cluster
          n <- sum(x_total$Freq) - m                             # cells NOT IN cluster
          k <- sum(x_test$Freq)                                      # cells from test condition
          x <- c(0:k)                                           # cells both IN cluster and from test condition
          
          # Use the dhyper built-in function for hypergeometric density 
          probabilities <- dhyper(x, m, n, k, log = FALSE)
          probs[clust,1:length(probabilities)] <- probabilities
        }
        
        probs <- as.data.frame(t(probs))
        probs$x <- x
        
        probs <- pivot_longer(probs, cols = 1:length(levels(x_test$Var1)), names_to = "cluster", values_to = "Prob")
        probs$cluster <- factor(probs$cluster)
        
        #Now calculate p-values from these probability distributions
        sig_test <- data.frame("unique" = rep(x_test$Var2[1], length(levels(probs$cluster))), "cluster" = levels(probs$cluster), "pval" = NA)
        rownames(sig_test) <- levels(probs$cluster)
        
        for (clust in sig_test$cluster) {
          #Cells from test condition and in cluster
          xIN <- x_test$Freq[which(x_test$Var1 == clust)]
          
          if (greater_than) {
            probs2 <- subset(probs, cluster == clust & x >= xIN)
          } else {
            probs2 <- subset(probs, cluster == clust & x <= xIN)
          }
          
          pval <- sum(probs2$Prob)
          
          sig_test[clust,"pval"] <- pval
          
        }
        
        return(sig_test)
        
      }, x_control = counts_split[[ctrl]])
      
      results <- bind_rows(results)
      return(results)
      
    }))
    
  } else {
    
    ctrl <- control_identifier
    
    #Create counts table of clustering assignments
    counts_true <- as.data.frame(table(df[,"cluster"], df[,group_cells_by]))
    counts_split <- split(counts_true, counts_true$Var2)
    
    #Iterate through each condition to get a p-value from Fisher's exact test for significant enrichment 
    #of cells from that condition in each cluster relative to the control condition
    results <- pblapply(counts_split, FUN = function(x_test, x_control) {
      
      #Initialize dataframe to store probabilities for each cluster
      probs <- data.frame()
      
      x_total <- bind_rows(x_test,x_control)
      
      #Now iterate over each cluster and get probabilities and p-values
      for (clust in levels(x_test$Var1)) {
        
        # Initialize variables
        m <- sum(x_total$Freq[which(x_total$Var1 == clust)])    # cells IN cluster
        n <- sum(x_total$Freq) - m                             # cells NOT IN cluster
        k <- sum(x_test$Freq)                                      # cells from test condition
        x <- c(0:k)                                           # cells both IN cluster and from test condition
        
        # Use the dhyper built-in function for hypergeometric density 
        probabilities <- dhyper(x, m, n, k, log = FALSE)
        probs[clust,1:length(probabilities)] <- probabilities
      }
      
      probs <- as.data.frame(t(probs))
      probs$x <- x
      
      probs <- pivot_longer(probs, cols = 1:length(levels(x_test$Var1)), names_to = "cluster", values_to = "Prob")
      probs$cluster <- factor(probs$cluster)
      
      #Now calculate p-values from these probability distributions
      sig_test <- data.frame("unique" = rep(x_test$Var2[1], length(levels(probs$cluster))), "cluster" = levels(probs$cluster), "pval" = NA)
      rownames(sig_test) <- levels(probs$cluster)
      
      for (clust in sig_test$cluster) {
        #Cells from test condition and in cluster
        xIN <- x_test$Freq[which(x_test$Var1 == clust)]
        
        if (greater_than) {
          probs2 <- subset(probs, cluster == clust & x >= xIN)
        } else {
          probs2 <- subset(probs, cluster == clust & x <= xIN)
        }
        
        pval <- sum(probs2$Prob)
        
        sig_test[clust,"pval"] <- pval
        
      }
      
      return(sig_test)
      
    }, x_control = counts_split[[ctrl]])
    
    results2 <- bind_rows(result)
    
  }
  
  results2$logP <- round((log10(results2$pval)*-1),6)
  results2 <- results2[,-3]
  results2 <- as.data.frame(pivot_wider(results2, names_from = cluster, values_from = logP))
  rownames(results2) <- results2[,1]
  results2 <- results2[,-1]
  results2 <- results2[,as.character(sort(as.numeric(colnames(results2))))]
  
  return(results2)
  
}



#Plotting function for Fisher exact test results
plot_cluster_enrichment <- function(test_results, 
                                    annotation_bar_groupings = NULL,
                                    annotation_colors = NULL,
                                    label_size = 6,
                                    cluster_column_var = TRUE,
                                    row_label_annotation_pval_thresh = 0.05) {
  
  #Determine which rows to annotate
  row_maximums <- rowMax(test_results)
  rows_to_label_pos <- which(row_maximums >= (log10(row_label_annotation_pval_thresh) * -1))
  rows_to_label <- rownames(test_results)[rows_to_label_pos]
  
  #Add annotation color bar
  if(!is.null(annotation_bar_groupings)) {
    
    classifier <- annotation_bar_groupings
    
    #Set color key if none provided
    if (is.null(annotation_colors)) {
      annotation_colors = colorRamp2(1:length(unique(classifier)), sample(gg_color_hue(length(unique(classifier)))))
    }
    
    row_annot <- rowAnnotation(Group = factor(classifier),
                               Cond = anno_mark(at= rows_to_label_pos, labels = rows_to_label, which="row", side = "right", labels_gp = gpar(col= "black", fontsize = label_size)),
                               col = list(Group = annotation_colors),
                               simple_anno_size = unit(2, "mm"),
                               show_annotation_name = F, show_legend = T)
    
    
    hm <- ComplexHeatmap::Heatmap(as.matrix(test_results),
                                  col = colorRamp2(c(0, 1.3, 1.4, 7), c("grey", "grey", "red", "darkred")),
                                  column_title = "cluster",
                                  column_title_side = "bottom",
                                  cluster_columns = cluster_column_var,
                                  row_title = NULL,
                                  #row_title_side = "right",
                                  show_row_names = F,
                                  row_dend_reorder = TRUE,
                                  column_names_gp = gpar(fontsize = label_size),
                                  #row_names_gp = gpar(fontsize = label_size),
                                  right_annotation = row_annot,
                                  name = "-log(pval)")
    
  } else {
    
    row_annot <- rowAnnotation(Cond = anno_mark(at= rows_to_label_pos, labels = rows_to_label, which="row", side = "right", labels_gp = gpar(col= "black", fontsize = label_size)),
                               show_annotation_name = F, show_legend = T)
    
    hm <- ComplexHeatmap::Heatmap(as.matrix(test_results),
                                  #col = colorRamp2(c(0, 1.3, 7), c("blue", "grey", "red")),
                                  col = colorRamp2(c(0, 1.3, 1.4, 7), c("grey", "grey", "red", "darkred")),
                                  column_title = "cluster",
                                  column_title_side = "bottom",
                                  cluster_columns = cluster_column_var,
                                  row_title = NULL,
                                  #row_title_side = "right",
                                  show_row_names = F,
                                  row_dend_reorder = TRUE,
                                  column_names_gp = gpar(fontsize = label_size),
                                  #row_names_gp = gpar(fontsize = label_size),
                                  right_annotation = row_annot,
                                  name = "-log(pval)")
  }
  
  return(hm)
  
}

#New plotting function for Fisher exact test results in new wide data format
plot_cluster_enrichment2 <- function(test_results, 
                                     annotation_colors,
                                     column_split_groups,
                                     col_map = NULL,
                                     label_size = 6) {
  
  n_clusters <- ncol(test_results)/length(annotation_colors)
  
  col_annot <- columnAnnotation(Base = rep(names(annotation_colors), each = n_clusters),
                                col = list(Base = annotation_colors),
                                simple_anno_size = unit(2, "mm"),
                                show_annotation_name = F, show_legend = T)
  
  if (is.null(col_map)) {
    col_map <- colorRamp2(c(0, 1.3, 2.7, 7), c("white", "grey", "red", "darkred"))
  }
  
  hm <- ComplexHeatmap::Heatmap(as.matrix(test_results),
                                col = col_map,
                                column_title = "cluster",
                                column_title_side = "bottom",
                                cluster_columns = F,
                                column_split = column_split_groups,
                                column_gap = unit(5, "mm"),
                                row_title = NULL,
                                show_row_names = T,
                                row_dend_reorder = TRUE,
                                column_names_gp = gpar(fontsize = label_size),
                                row_names_gp = gpar(fontsize = label_size),
                                top_annotation = col_annot,
                                border = T,
                                name = "-log10(pval)")
  
  return(hm)
  
}



## Plot single cell scores for genesets
##   df: a data frame.
##   set: column name of the geneset score to plot
plot_geneset_scores <- function(df, set, label_size = 6, plot_type = c("ridgeline", "boxplot", "UMAP")) {
  
  df$median <- ave(df[,set], as.factor(df[,"unique"]), FUN=median)
  
  if (plot_type == "ridgeline") {
    
    p <- ggplot(df, aes_string(x = set, y = "unique", fill = "median")) +
      geom_density_ridges(rel_min_height = 0.01) +
      scale_fill_viridis(option = "C", direction = 1) +
      xlim(quantile(geneset_scores[,set], probs = 0.01), quantile(geneset_scores[,set], probs = 0.99)) +
      facet_wrap(~base_new, scales = "free_y") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = label_size))
    
  } else if (plot_type == "boxplot") {
    
    p <- ggplot(df, aes_string(x = set, y = "unique", color = "median")) +
      geom_boxplot() +
      scale_color_viridis(option = "C", direction = 1) +
      facet_wrap(~base_new, scales = "free_y") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = label_size))
    
  } else if (plot_type == "UMAP") {
    
    p <- ggplot(df, aes_string(x = "umap_1", y = "umap_2", color = set)) +
      geom_point(size = 0.3, alpha = 1) +
      scale_color_viridis(option = "C", direction = 1) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
  }
  
  return(p)
  
}

## Test for differential expression of genesets
##   df: a data frame of geneset scores.
##   test: which statistical test to run, KS or Wilcoxin
##   control: string indicating control condition to test against

##MJW - update to allow for splitting by base condition. Take method from Fisher test function
geneset_DE_test <- function(df, test = c("KS", "Wilcoxin"), control) {
  
  #Group cells into conditions
  df_split <- split(df, df$unique)
  
  #control arg must match one of the names of df_split
  assertthat::assert_that(control %in% names(df_split),
                          msg = paste('Control condition must be one of the conditions in the dataframe'))
  
  #Parse the columns to keep
  n_sets <- length(colnames(df)) - 4
  control_group <- df_split[[control]]
  
  result <- lapply(df_split, FUN = function(x) {
    
    res_out <- mapply(FUN = function(df_test,df_ctrl) {
      if (test == "KS") { 
        res <- suppressWarnings(ks.test(df_test, df_ctrl)) 
      } else { 
        res <- wilcox.test(df_test, df_ctrl, paired = F)
      }
      
      return(res$p.value)
    }, df_test = x[,1:n_sets], df_ctrl = control_group[,1:n_sets])
    
    return(res_out)
    
  })
  
  test_results <- as.matrix(do.call(cbind, result))
  test_results <- log10(test_results) * -1
  return(test_results)
  
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



#LASSO Regression functions (Cao et al.)
do_lasso <- function(gene_matrix, gene_vector, seed = 1) {
  
  M1 = gene_matrix
  y = gene_vector
  
  set.seed(seed)
  cv.out1 = cv.glmnet(M1, y, alpha=1, lambda=exp(seq(log(0.001), log(10), length.out=100)))
  r2_1 = r2_glmnet(cv.out1, y)
  
  bestlam = cv.out1$lambda.1se
  cor_list = coef(cv.out1, s= bestlam)
  cor_length = length(cor_list)
  df_cor = data.frame("id" = row.names(cor_list)[2:cor_length], "corcoef" = cor_list[2:cor_length])
  return(list(r2_1, df_cor))
}

r2_glmnet <- function(cv.out, y) {
  bestlam = cv.out$lambda.1se
  i = which(cv.out$lambda == bestlam)
  e <- cv.out$cvm[i]
  r2 <- 1 - e / var(y)
  if(r2 < 0)
  {
    r2 = 0
  }
  return(r2)
}

#Matt's
create_lasso_ggc_matrix <- function(lasso_output, gene_list, symmetrize = TRUE) {
  
  # #Drop coefficients into a new correlation matrix
  # ggc_mat <- matrix(nrow = length(gene_list), ncol = length(gene_list))
  # rownames(ggc_mat) <- gene_list
  # colnames(ggc_mat) <- gene_list
  # 
  # for (i in 1:length(lasso_output)) {
  #   result <- lasso_output[[i]][[2]]
  #   row_match <- match(result$id, rownames(ggc_mat))
  #   ggc_mat[row_match,i] <- result$corcoef
  # }
  
  ggc_mat <- pblapply(1:length(lasso_output), function(x) {
    result <- lasso_output[[x]][[2]]
    
    #Add back in self correlation
    result <- rbind(c(gene_list[x],NA), result)
    
    #Put into order
    result$index <- match(result$id, gene_list)
    result <- result[order(result$index, decreasing = F),]
    
    return(as.numeric(result$corcoef))
    
  })
  
  #Create matrix
  ggc_mat <- as.matrix(do.call(cbind, ggc_mat))
  rownames(ggc_mat) <- colnames(ggc_mat) <- gene_list
  
  #Filter genes with no correlation
  coeff_sums <- Matrix::colSums(abs(ggc_mat), na.rm = T)
  top_genes <- rownames(ggc_mat)[which(coeff_sums > 0.00)]
  
  ggc_mat <- ggc_mat[top_genes,top_genes]
  
  if(symmetrize) {
    ggc_mat = (ggc_mat+t(ggc_mat))/2   # symmetrize the coeff matrix
  }
  
  return(ggc_mat)
  
}

bootstrap_ggc_matrix <- function(expr_matrix, sampling_iterations = 10) {
  
  all_genes <- colnames(expr_matrix)
  barcodes <- rownames(expr_matrix)
  
  #Sample 70% of the cells in the expression matrix containing the ~1800 genes in the lasso matrix
  #Make pearson correlation matrix of the genes
  message("Sampling Data")
  sampled_matrices <- pbmclapply(1:sampling_iterations, function(i) {
    samp <- sample(barcodes, (0.7*length(barcodes)), replace = F)
    new_mat <- expr_matrix[samp,]
    cor_mat <- cor(new_mat)
    return(cor_mat)
  }, mc.cores = 2)
  
  return(sampled_matrices)
  
}

bootstrap_apclust <- function(ggc_matrices, q_range) {
  
  ap.out <- lapply(q_range, FUN = function(q) { 
    
    message("Peforming clustering for q = ", q)
    output <- pbmclapply(ggc_matrices, function(mat) {
      res <- apcluster::apcluster(s = mat, q = q)
      return(res)
    }, mc.cores = 2)
    
    return(output)
  })
  
  return(ap.out)
  
}

make_consensus_matrix <- function(cluster_results, ggc_matrix) {
  
  #Initialize zero matrices for counts
  counts.paired <- ggc_matrix
  counts.paired[,] <- 0
  counts.total <- counts.paired
  
  message("Generating consensus matrix")
  for (res in cluster_results) {
    res.clust <- res@clusters
    #Get all genes sampled in this sub_matrix
    all_sampled <- names(unlist(res.clust))
    
    #Count all pairs sampled together
    counts.total[all_sampled,all_sampled] <- counts.total[all_sampled,all_sampled] + 1
    
    #Count all pairs clustered together
    for (c in res.clust) {
      ids <- names(c)
      counts.paired[ids,ids] <- counts.paired[ids,ids] + 1
    }
  }
  
  consensus.out <- counts.paired/counts.total
  return(consensus.out)
  
}

plot_consensus_matrix <- function(mat) {
  d <- distanceMatrix(mat, metric = 'spearman')
  d <- hclust(d, method = "average")
  set.seed(12)
  d <- reorder(as.dendrogram(d), rowMeans(mat, na.rm = T))
  mat2 <- mat[labels(d),labels(d)]
  hm <- Heatmap(mat2,
                show_row_dend = F,
                show_column_dend = F,
                show_column_names = F,
                show_row_names = F,
                cluster_columns = F,
                cluster_rows = F)
  return(hm)
}

evaluate_clustering <- function(sim_matrix, ggc_mat, method = c("hclust", "AP"), min_size = 3) {
  
  if (method == "AP") {
    
    ap.out <- apcluster::apcluster(s = sim_matrix, q = 0)
    clusters <- ap.out@clusters
    #Create data frame to store module data
    modules <- bind_rows(mapply(FUN = function(x, i) {
      out_df <- data.frame(module = i, gene = names(x))
      return(out_df)
      
    }, x = clusters, i = 1:length(clusters), SIMPLIFY = F))
    
  } else {
    
    d <- hclust(dist(sim_matrix))
    h_cutoff <- quantile(d$height, 0.99)
    
    modules <- data.frame("module" = cutree(d, h = h_cutoff, order_clusters_as_data = F))
    modules$gene <- rownames(modules)
  }
  
  #Generate random clusters
  modules$gene_random <- sample(modules$gene)
  
  mods <- split(modules, modules$module)
  mod_size <- data.frame(module = 1:length(mods), n = sapply(mods, function(x) {length(x$gene)}))
  mods <- mods[mod_size$module[which(mod_size$n >= min_size)]]
  
  #Rename modules in case some were removed
  mods2 <- bind_rows(mapply(FUN = function(mod, name) {
    mod$module <- name
    return(mod)
  }, mod = mods, name = 1:length(mods), SIMPLIFY = F))
  
  colnames(mods2)[2:3] <- c("true", "random")
  
  mods2 <- pivot_longer(mods2, cols = 2:3, names_to = "set", values_to = "gene")
  
  DB <- lapply(split(mods2, mods2$set), function(df) {
    #Subset matrix on genes included in clusters
    cor_mat2 <- ggc_mat[df$gene,df$gene]
    indices <- match(df$gene, rownames(cor_mat2))
    df$index <- indices
    cl_df <- df[order(df$index, decreasing = F),]
    
    testDB <- index.DB(x = cor_mat2, cl = cl_df$module)
  })
  
  return(DB)
  
}

tune_clusters <- function(optimal_consensus_matrix, min_module_size = 5, min_n_clusters = 30, DB_prc = 0.99) {
  
  #Get initial clusters
  apclust <- apcluster::apcluster(s = optimal_consensus_matrix, q = 0)
  apclust2 <- apcluster::aggExCluster(s = optimal_consensus_matrix, x = apclust, includeSim = T)
  
  og_mods <- apclust@clusters
  mods <- og_mods
  minimum <- min(sapply(mods, function(x) {length(x)}))
  
  #Merge clusters until the smallest module is above the minimum size threshold while maintaining a minimum number of total clusters
  
  while (minimum < min_module_size) {
    current_k <- length(mods)
    new_k <- current_k-1
    merge <- apcluster::cutree(apclust2, k=new_k)
    mods <- merge@clusters
    minimum <- min(sapply(mods, function(x) {length(x)}))
    
    #If this minimum resulted in too few clusters, reduce the minimum module size and restart
    if (minimum >= min_module_size & length(mods) < min_n_clusters) {
      mods <- og_mods
      minimum <- min(sapply(mods, function(x) {length(x)}))
      min_module_size <- min_module_size-1
    }
  }
  
  modules <- bind_rows(mapply(FUN = function(x, i) {
    out_df <- data.frame(module = i, gene_index = x)
    return(out_df)
    
  }, x = mods, i = 1:length(mods), SIMPLIFY = F))
  
  modules$gene <- rownames(optimal_consensus_matrix)[modules$gene_index]
  
  #DB index
  cor_mat2 <- ggc_mat[modules$gene,modules$gene]
  indices <- match(modules$gene, rownames(cor_mat2))
  modules$index2 <- indices
  cl_df <- modules[order(modules$index2, decreasing = F),]
  
  DB_start <- index.DB(x = cor_mat2, cl = cl_df$module)
  DB_start <- DB_start$DB
  
  #Now remove individual clusters and see if the DB index improves
  modules2 <- modules
  
  DB_old <- DB_start
  DB_new <- DB_old*DB_prc
  while(DB_new <= DB_old*DB_prc) {
    DB_old <- DB_new
    
    DBtest <- pbsapply(1:length(unique(modules2$module)), function(x) {
      modules3 <- subset(modules2, module != x)
      mods3_split <- split(modules3, modules3$module)
      modules3 <- bind_rows(mapply(FUN = function(mod, name) {
        mod$module <- name
        return(mod)
      }, mod = mods3_split, name = 1:length(mods3_split), SIMPLIFY = F))
      
      cor_mat2 <- ggc_mat[modules3$gene,modules3$gene]
      indices <- match(modules3$gene, rownames(cor_mat2))
      modules3$index2 <- indices
      cl_df <- modules3[order(modules3$index2, decreasing = F),]
      DB2 <- index.DB(x = cor_mat2, cl = cl_df$module)
      return(DB2$DB)
    })
    
    if(any(DBtest < (DB_old*DB_prc))) {
      DB_new <- min(DBtest)
      mod_to_remove <- which(DBtest == DB_new)
      modules2 <- subset(modules2, module != mod_to_remove)
      mods2_split <- split(modules2, modules2$module)
      modules2 <- bind_rows(mapply(FUN = function(mod, name) {
        mod$module <- name
        return(mod)
      }, mod = mods2_split, name = 1:length(mods2_split), SIMPLIFY = F))
    }
    
  }
  
  return(modules2)
  
}


#Set ggplot theme
theme_set(theme_bw())


# # Matt's functions for machine learning classification models
# #Probably need to be updated for this project
# make_kfold_splits <- function(df, kfolds = 5, ntimes = 5) {
#   
#   #Folds are the number of equal-sized partitions of the dataframe after stratification,
#   #preserving relative frequencies of each condition per fold.
#   
#   #Times is how many times to independently repeat the k-fold splitting
#   
#   splits <- caret::createMultiFolds(df$Cond, k = kfolds, times = ntimes)
#   splits
#   
# }
# 
# do_RF <- function(ksplit, df, class_name = "Cond", predictor_vars, mtry = NA, tune_mtry = FALSE, downsample = T) {
#   
#   #Split data
#   training_data <- df[ksplit,]
#   test_data <- df[-ksplit,]
#   
#   #Train model
#   if (is.na(mtry)) {
#     
#     if(tune_mtry) {
#       #Optimize mtry on error rate
#       rf <- tuneRF(training_data[,predictor_vars], training_data[,class_name],
#                    stepFactor = 0.6,
#                    plot = FALSE,
#                    ntreeTry = 150,
#                    trace = TRUE,
#                    improve = 0.05,
#                    doBest = TRUE)
#     } else {
#       
#       #Write formula
#       rf_formula <- as.formula(paste(paste0(class_name, " ~"), paste(predictor_vars, collapse = " + ")))
#       #Use default mtry value
#       rf <- randomForest(rf_formula, data = training_data, proximity = TRUE, na.action = na.roughfix)
#       
#     }
#     
#     
#   } else {
#     
#     #Use specified mtry value
#     rf_formula <- as.formula(paste(paste0(class_name, " ~"), paste(predictor_vars, collapse = " + ")))
#     rf <- randomForest(rf_formula, data = training_data, mtry = mtry, proximity = TRUE, na.action = na.roughfix)
#     
#   }
#   
#   #Validate model
#   #Downsample validation data for balanced test data
#   if (downsample) {
#     test_data <- caret::downSample(test_data[, colnames(test_data)[!(colnames(test_data) == "Cond")] ], test_data$Cond, yname = "Cond")
#   }
#   
#   p1 <- predict(rf, test_data)
#   p1_output <- confusionMatrix(p1, test_data[,class_name], mode = "everything")
#   p1_output
#   
# }
# 
# get_confusion_matrix <- function(result, apply_dim) {
#   conf_mat <- result$table
#   CM <- matrix(conf_mat, ncol = ncol(conf_mat), dimnames = dimnames(conf_mat))
#   
#   if (apply_dim == 1) {
#     CM <- matrix(apply(CM, apply_dim, compute_percentage), ncol = ncol(conf_mat), nrow = ncol(conf_mat), dimnames = dimnames(CM), byrow = T)
#   } else if (apply_dim == 2) {
#     CM <- matrix(apply(CM, apply_dim, compute_percentage), ncol = ncol(conf_mat), nrow = ncol(conf_mat), dimnames = dimnames(CM), byrow = F)
#   }
#   
#   CM
# }
# 
# get_model_performance <- function(result) {
#   as.matrix(result$byClass)
# }
# 
# get_model_stats <- function(result) {
#   result$overall
# }
# 
# summarize_kfold_model <- function(model_results, apply_dim) {
#   
#   summary <- list()
#   
#   l <- length(model_results)
#   
#   #Average all confusion matrices
#   CM_list <- lapply(model_results, FUN = get_confusion_matrix, apply_dim = apply_dim)
#   CM_mean <- Reduce("+", CM_list)/l
#   summary[["CM"]] <- CM_mean
#   
#   #mean +/- SD for performance metrics
#   metrics_list <- lapply(model_results, FUN = get_model_performance)
#   
#   metrics_mean <- as.data.frame(Reduce("+", metrics_list)/l)
#   metrics_sd <- apply(array(unlist(metrics_list), c(dim(metrics_mean), l)), c(1,2), sd)
#   
#   metrics_mean$Cond <- factor(gsub("Class: ", '', rownames(metrics_mean)), levels = rownames(CM_mean))
#   metrics_mean <- pivot_longer(metrics_mean, cols = -(length(colnames(metrics_mean))), names_to = "metric", values_to = "score")
#   
#   metrics_sd <- as.data.frame(metrics_sd)
#   colnames(metrics_sd) <- metrics_mean$metric[1:11]
#   metrics_sd$Cond <- levels(metrics_mean$Cond)
#   metrics_sd <- pivot_longer(metrics_sd, cols = 1:11, names_to = "metric", values_to = "SD")
#   
#   summary[["metrics"]] <- as.data.frame(metrics_mean)
#   summary[["metrics_sd"]] <- as.data.frame(metrics_sd)
#   
#   #overall statistics
#   overall_stats <- lapply(model_results, FUN = get_model_stats)
#   overall_mean <- Reduce("+", overall_stats)/l
#   summary[["overall"]] <- overall_mean
#   
#   summary
#   
# }
# 
# make_CM_heatmap <- function(mod, breaks = c(0,40,80)) {
#   labels_rows <- rowAnnotation(Cond = factor(rownames(mod$CM), levels = rownames(mod$CM)),
#                                col = list(Cond = pMHC_colors),
#                                simple_anno_size = unit(2, "mm"),
#                                show_annotation_name = FALSE, show_legend = F)
#   
#   labels_columns <- HeatmapAnnotation(Cond = factor(rownames(mod$CM), levels = rownames(mod$CM)),
#                                       col = list(Cond = pMHC_colors),
#                                       simple_anno_size = unit(2, "mm"),
#                                       show_annotation_name = FALSE, show_legend = F)
#   
#   CM_palette <- colorRamp2(breaks, c("blue", "white", "red"))
#   
#   hm <- Heatmap(mod$CM, cluster_rows = F, cluster_columns = F,
#                 col <- CM_palette,
#                 cell_fun = function(j, i, x, y, width, height, fill) {
#                   grid.text(paste0(round(mod$CM[i, j],1), "%"), x, y, gp = gpar(fontsize = 10, fontface = "bold"))},
#                 name = "Pred. %",
#                 row_title = "Predicted",
#                 right_annotation = labels_rows,
#                 bottom_annotation = labels_columns,
#                 show_row_names = F, show_column_names = F,
#                 column_title = "True")
#   
#   return(hm)
#   
# }
# 
# make_CM_heatmap2 <- function(mod, breaks = c(0,40,80)) {
#   
#   CM_palette <- colorRamp2(breaks, c("blue", "white", "red"))
#   
#   hm <- Heatmap(mod$CM, cluster_rows = F, cluster_columns = F,
#                 col <- CM_palette,
#                 name = "Pred. %",
#                 row_title = "Predicted",
#                 show_row_names = T, show_column_names = T,
#                 column_title = "True")
#   
#   return(hm)
#   
# }
# 
# get_stat_byClass <- function(m, name, stat = "F1") {
#   out.df <- subset(m$metrics, metric == stat)
#   sd <- subset(m$metrics_sd, metric == stat)
#   out.df$SD <- sd$SD
#   out.df$model <- as.character(name)
#   as.data.frame(out.df)
# }
# 
# 
# #Other functions
# optimize_lasso_APclustering <- function(results, expr_mat, optimization_metric = "cov", random_control = TRUE, minimum_module_size = 2) {
#   
#   mods <- results@clusters
#   
#   #Remove clusters with fewer than a threshold number of genes
#   mod_size <- sapply(mods, function(x) {length(x)})
#   mods <- mods[which(mod_size >= minimum_module_size)]
#   
#   #Extract gene lists
#   gene_lists <- lapply(mods, function(x) {names(x) })
#   
#   #Now subset expression matrix on each set of genes and compute optimization metric
#   test_metrics <- lapply(gene_lists, function(x) {
#     mat <- expr_mat[,x]
#     
#     if (optimization_metric == "cov") {
#       output <- mean(stats::cov(mat))
#     } else {
#       output <- mean(apply(mat, 1, optimization_metric, na.rm=TRUE))
#     }
#     
#     return(output)
#   })
#   
#   
#   if (random_control) {
#     
#     module_df <- bind_rows(mapply(FUN = function(m, i) {
#       out_df <- data.frame(index = m, module = i, gene = names(m))
#       return(out_df)
#       
#     }, m = mods, i = 1:length(mods), SIMPLIFY = F))
#     
#     module_df$gene_random <- sample(module_df$gene)
#     gene_lists_random <- lapply(split(module_df, module_df$module), function(x) { x$gene_random })
#     
#     test_random <- lapply(gene_lists_random, function(y) {
#       mat <- expr_mat[,y]
#       
#       if (optimization_metric == "cov") {
#         output <- mean(stats::cov(mat))
#       } else {
#         output <- mean(apply(mat, 1, optimization_metric, na.rm=TRUE))
#       }
#       
#       return(output)
#     })
#     
#   }
#   
#   out_df <- data.frame( mean = c(mean(unlist(test_metrics)), mean(unlist(test_random))),
#                         var = c(stats::var(unlist(test_metrics)), stats::var(unlist(test_random))),
#                         modules = c("test", "random") )
#   
#   return(out_df)
#   
# }
# 
# Tune_APclustering <- function(k_range, similarity_matrix, sc_expression_matrix, minimum_module_size = 2) {
#   
#   #Perform clustering with variable number of clusters
#   print("Performing clustering")
#   ap_test <- pbmclapply(k_range, FUN = function(k) { apcluster::apclusterK(s = similarity_matrix, K = k, prc = 5, bimaxit = 20, verbose = T) }, mc.cores = 2)
#   
#   #Compute degree of co-expression within each cluster
#   print("Computing intra-cluster covariances")
#   ap_opt <- pbmclapply(ap_test, FUN = optimize_lasso_APclustering, expr_mat = sc_expression_matrix, mc.cores = 2)
#   
#   #Output data frame for plotting
#   ap_opt <- bind_rows(ap_opt)
#   ap_opt$k <- rep(unlist(lapply(ap_test, function(x) { length(x@clusters) })), each = 2)
#   ap_opt <- bind_rows(lapply(split(ap_opt, ap_opt$modules), function(df) { df[match(unique(df$k), df$k),] }))
#   
#   return(list(ap_test, ap_opt))
# }
# 
# export_heatmap <- function(hm, name, out_dir) {
#   setwd(out_dir)
#   tiff(file = paste0(as.character(name),".tiff"), width = 360, height = 360)
#   draw(hm)
#   dev.off()
# }
# 
# tune_sc_groups <- function(agg_mat, apclust_result, DB_prc = 0.99) {
#   
#   #Get initial clusters
#   apclust <- apclust_result[[2]]
#   apclust2 <- apclust_result[[3]]
#   
#   og_groups <- apclust@clusters
#   groups <- og_groups
#   
#   merged_groups <- bind_rows(mapply(FUN = function(x, i) {
#     out_df <- data.frame(group = i, sc_index = x)
#     return(out_df)
#     
#   }, x = groups, i = 1:length(groups), SIMPLIFY = F))
#   
#   #Save starting groups
#   merged_groups2 <- merged_groups
#   merged_groups2$barcode <- rownames(merged_groups2)
#   
#   #Reduce the number of clusters by 1 and test DB index
#   agg_mat2 <- agg_mat[,merged_groups2$barcode]
#   
#   DB_start <- index.DB(x = t(agg_mat2), cl = merged_groups2$group)
#   DB_start <- DB_start$DB
#   
#   DB_old <- DB_start
#   DB_new <- DB_old*DB_prc
#   
#   while(DB_new <= DB_old*DB_prc) {
#     DB_old <- DB_new
#     current_k <- length(groups)
#     new_k <- current_k-1
#     merge <- apcluster::cutree(apclust2, k=new_k)
#     groups <- merge@clusters
#     merged_groups2 <- bind_rows(mapply(FUN = function(x, i) {
#       out_df <- data.frame(group = i, sc_index = x)
#       return(out_df)
#       
#     }, x = groups, i = 1:length(groups), SIMPLIFY = F))
#     
#     #Cutree removes some of the row naems so need to update them from the original matrix
#     merged_groups2$barcode <- colnames(agg_mat)[merged_groups2$sc_index]
#     
#     agg_mat2 <- agg_mat[,merged_groups2$barcode]
#     
#     DB2 <- index.DB(x = t(agg_mat2), cl = merged_groups2$group)
#     DB_new <- DB2$DB
#     
#   }
#   
#   return(merged_groups2)
#   
# }
# 
# tune_sc_groups2 <- function(agg_mat, apclust_result, k) {
#   
#   #Get initial clusters
#   apclust <- apclust_result[[2]]
#   apclust2 <- apclust_result[[3]]
#   
#   merge <- apcluster::cutree(apclust2, k=k)
#   groups <- merge@clusters
#   merged_groups2 <- bind_rows(mapply(FUN = function(x, i) {
#     out_df <- data.frame(group = i, sc_index = x)
#     return(out_df)
#     
#   }, x = groups, i = 1:length(groups), SIMPLIFY = F))
#   
#   #Cutree removes some of the row naems so need to update them from the original matrix
#   merged_groups2$barcode <- colnames(agg_mat)[merged_groups2$sc_index]
#   
#   return(merged_groups2)
#   
# }



