#Start with the cds containing CD3/28 base condition and only highest dose of each cytokine
#Aggregate by gene module and condition - outputting z-scores
agg <- sc_pseudobulk(cds_hi_dose[[3]], gene_groups = output_modules[,c("gene", "module")], zscore = T, cell_group_column = "Cond")

s <- negDistMat(t(agg))
#Run clustering
ap_cyt <- apcluster(s = s, q = 0.6)
groups <- ap_cyt@clusters
ap_cyt2 <- aggExCluster(s = s, x = ap_cyt, includeSim = T)
apcluster::heatmap(ap_cyt2)

cyt_order <- names(unlist(groups[ap_cyt2@order]))

#Create data frame to store cytokine groups
cyt_groups <- bind_rows(mapply(FUN = function(x, i) {
  out_df <- data.frame(index = x, group = i, cytokine = names(x))
  return(out_df)
  
}, x = groups, i = 1:length(groups), SIMPLIFY = F))

#reorder matrix columns
agg_norm3 <- agg[,cyt_order]

#Add bar annotations to identify cytokine groups
#Sort data frame to match matrix order
cyt_groups_ordered <- cyt_groups[cyt_order,]
cyt_groups_ha <- HeatmapAnnotation(Module = cyt_groups_ordered$group,
                                   col = list(Module = colorRamp2(1:length(groups), gg_color_hue(length(groups)))),
                                   simple_anno_size = unit(2, "mm"),
                                   show_annotation_name = FALSE, show_legend = F)

col_split <- data.frame(split = factor(cyt_groups_ordered$group))

Heatmap(agg_norm3,
              #col = colorRamp2(c(-1,-0.5,0,0.5,1), c("darkblue", "blue", "white", "red", "darkred")),
              show_row_dend = F,
              show_column_dend = T,
              show_column_names = T,
              column_names_gp = gpar(fontsize = 6),
              row_names_gp = gpar(fontsize = 6),
              show_row_names = T,
              row_names_side = "left",
              cluster_columns = T,
              cluster_rows = T,
              bottom_annotation = cyt_groups_ha,
              #right_annotation = rowAnnotation(FCrange = anno_barplot(gm_range), annotation_name_rot = 90, gap = unit(2, "mm")),
              column_split = col_split,
              column_title_side = "bottom",
              border = T,
              row_title = "Gene Module",
              name = "z-score")


#################################################################

d <- distanceMatrix(agg, metric = 'spearman')
d <- hclust(d, method = "average")
set.seed(12)
d <- reorder(as.dendrogram(d), rowMeans(t(agg), na.rm = T))

groups <- cutree(d, k = 8)
col_split <- data.frame(split = factor(groups))


Heatmap(agg,
        show_row_dend = F,
        show_column_dend = T,
        show_column_names = T,
        column_names_gp = gpar(fontsize = 5),
        row_names_gp = gpar(fontsize = 8),
        show_row_names = T,
        row_names_side = "left",
        cluster_columns = d,
        cluster_rows = T,
        #bottom_annotation = cyt_groups_ha,
        #right_annotation = rowAnnotation(FCrange = anno_barplot(gm_range), annotation_name_rot = 90, gap = unit(2, "mm")),
        column_split = 10,
        column_title_side = "bottom",
        border = T,
        row_title = "Gene Module",
        name = "z-score")


