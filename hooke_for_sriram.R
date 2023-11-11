#' The cell_count_set class
#'
#' The main class used by Hooke to hold cell abundances data.
#' cell_count_set extends the Monocle's cell_data_set class.
#'
#' This class is initialized from a matrix of expression values along with cell
#' and feature metadata.
#'
#' @field cds cell_data_set, the Monocle cell data set object that this class models.
#' @name cell_count_set
#' @rdname cell_count_set
#' @aliases cell_count_set-class
#' @exportClass cell_count_set
setClass("cell_count_set",
         contains = "cell_data_set",
         slots = c(cds = "cell_data_set",
                   cds_coldata = "tbl_df",
                   cds_reduced_dims = "SimpleList",
                   info = "SimpleList")
)
setMethod("is.na", "cell_count_set", function(x) FALSE)

new_cell_count_set <- function(cds,
                               sample_group,
                               cell_group,
                               sample_metadata = NULL,
                               cell_metadata = NULL,
                               lower_threshold = NULL,
                               upper_threshold = NULL,
                               keep_cds=TRUE,
                               norm_method = c("size_factors","TSS", "CSS",
                                               "RLE", "GMPR", "Wrench", "none"),
                               size_factors = NULL,
                               pseudocount = 0) {
  
  assertthat::assert_that(is(cds, 'cell_data_set'),
                          msg = paste('Argument cds must be a cell_data_set.'))
  
  assertthat::assert_that(sample_group %in% colnames(colData(cds)),
                          msg = paste('Argument sample_group value must be a column name in the cell_data_set.'))
  
  assertthat::assert_that(cell_group %in% colnames(colData(cds)),
                          msg = paste('Argument cell_group value must be a column name in the cell_data_set.'))
  
  assertthat::assert_that(is.null(sample_metadata) || is.data.frame(sample_metadata),
                          msg = paste('Argument sample_metadata must be a data frame.'))
  
  sample_group_names_cds <- unique(colData(cds)[[sample_group]])
  assertthat::assert_that(is.null(sample_metadata) || nrow(sample_metadata) == length(sample_group_names_cds),
                          msg = paste('Argument sample_metadata must have the same',
                                      'number of rows as there are distinct sample',
                                      'names in the cds column data.'))
  
  assertthat::assert_that(is.null(sample_metadata) || all(sample_group_names_cds %in% sample_metadata[['sample']]),
                          msg = paste('Argument sample_metadata must have sample group names in',
                                      'a column called \'sample\'.'))
  
  assertthat::assert_that(is.null(cell_metadata) || is.data.frame(cell_metadata),
                          msg = paste('Argument cell_metadata must be a data frame.'))
  
  cell_group_names_cds <- unique(colData(cds)[[cell_group]])
  assertthat::assert_that(is.null(cell_metadata) || nrow(cell_metadata) == length(cell_group_names_cds),
                          msg = paste('Argument cell_metadata must have the same',
                                      'number of rows as there are distinct cell_group',
                                      'names in the cds column data.'))
  
  assertthat::assert_that(is.null(cell_metadata) || all(cell_group_names_cds %in% row.names(cell_metadata)),
                          msg = paste('Argument cell_metadata row names must contain the cell_group',
                                      'names.'))
  
  assertthat::assert_that(is.null(lower_threshold) || is.numeric(lower_threshold),
                          msg = paste('Argument lower_threshold must be numeric.'))
  
  assertthat::assert_that(is.null(upper_threshold) || is.numeric(upper_threshold),
                          msg = paste('Argument upper_threshold must be numeric.'))
  
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(norm_method) == "", TRUE, TRUE),
             error = function(e) FALSE),
    msg = paste('Argument norm_method must be one of "size_factors",',
                '"TSS", "CSS", "RLE", "GMPR", "Wrench", or "none".'))
  norm_method <- match.arg(norm_method)
  
  if(sample_group != 'sample')
    colData(cds)$sample = NULL
  
  # check if anything contains NAs in it
  # if so drop them
  num_sample_group_NAs = sum(is.na(colData(cds)[[sample_group]]))
  if (num_sample_group_NAs != 0) {
    message(paste(num_sample_group_NAs, "NAs found in sample group. Dropping NAs."))
    cds = cds[, !is.na(colData(cds)[[sample_group]])]
  }
  
  num_cell_group_NAs = sum(is.na(colData(cds)[[cell_group]]))
  if (num_cell_group_NAs != 0) {
    message(paste(num_cell_group_NAs, "NAs found in cell group. Dropping NAs."))
    cds = cds[, !is.na(colData(cds)[[cell_group]])]
  }
  
  if (is.character(colData(cds)[[cell_group]])){
    num_cell_group_blanks = sum(colData(cds)[[cell_group]] == "")
    if (num_cell_group_blanks != 0) {
      message(paste(num_cell_group_blanks, "unlabeled cells found in cell group. Dropping unlabled cells."))
      cds = cds[, colData(cds)[[cell_group]] != ""]
    }
  }
  
  coldata_df = colData(cds) %>% tibble::as_tibble()
  # current commented out bc mess w projection clusters
  # coldata_df$cluster = monocle3::clusters(cds)
  # coldata_df$partition = partitions(cds)
  
  coldata_df = coldata_df %>% dplyr::rename("sample" = sample_group, "cell_group" = as.character(cell_group))
  #coldata_df$cell_group = factor(coldata_df$cell_group, levels=unique(colData(cds)[,cell_group]))
  
  coldata_df$group_id = coldata_df %>%
    dplyr::group_by(sample, cell_group) %>%
    dplyr::group_indices() %>% as.character
  
  # add to cds
  colData(cds)$group_id = coldata_df$group_id
  
  cds_summary = coldata_df %>%
    dplyr::group_by(sample, cell_group) %>%
    dplyr::summarize(cells = dplyr::n())
  
  cds_covariates_df = coldata_df %>%
    dplyr::select(-cell_group) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarize(across(where(is.numeric), function(x){mean(x)}),
                     across(where(is.factor), function(x) { tail(names(sort(table(x))), 1) }),
                     across(where(is.character), function(x) { tail(names(sort(table(x, useNA="ifany"))), 1) } ))
  
  if (is.null(sample_metadata) == FALSE){
    cds_covariates_df = left_join(cds_covariates_df, sample_metadata, by=c("sample"="sample"))
  }
  
  cds_covariates_df = cds_covariates_df %>% as.data.frame(cds_covariates_df, stringsAsFactors=FALSE)
  row.names(cds_covariates_df) = cds_covariates_df %>% dplyr::pull(sample)
  
  cell_counts_wide = tidyr::spread(cds_summary, sample, cells, fill=0)
  cell_states = as.character(cell_counts_wide %>% dplyr::pull("cell_group"))
  cell_counts_wide = as.matrix(cell_counts_wide[,2:ncol(cell_counts_wide)])
  
  row.names(cell_counts_wide) = cell_states
  
  # filter out cell groups based on counts
  if (is.null(lower_threshold) == FALSE) {
    cell_counts_wide = cell_counts_wide[Matrix::rowSums(cell_counts_wide) >= lower_threshold, ]
  }
  if (is.null(upper_threshold) == FALSE) {
    cell_counts_wide = cell_counts_wide[Matrix::rowSums(cell_counts_wide) <= upper_threshold, ]
  }
  
  # remove from cds
  removed_cell_states = setdiff(cell_states, rownames(cell_counts_wide))
  
  #cell_counts_wide = t(cell_counts_wide)
  
  cds_covariates_df = cds_covariates_df[colnames(cell_counts_wide),]
  
  # This is super confusing because of the way the arguments are
  # named in new_cell_data_set. We are making a matrix of
  # dimension MxN, where M are cell types and N are samples
  # (e.g. embryos, replicates, etc). The "gene" metadata monocle
  # normally expects will actually be used to hold cell group
  # metadata.
  
  ccs_cds = cds[, !colData(cds)[[cell_group]] %in% removed_cell_states]
  
  # TODO: We could probably avoid duplicating this info if keep_cds == TRUE, providing it
  # through accessor functions directly from the cds
  
  # FIXME: potentially we should be using the filtered one above? Potentially rename cell_group, sample, etc?
  cds_coldata = colData(ccs_cds) %>% as_tibble
  cds_reducedDims = reducedDims(ccs_cds)
  
  cell_metadata_subset <- cell_metadata[rownames(cell_counts_wide),,drop=FALSE]
  ccs = methods::new("cell_count_set",
                     monocle3::new_cell_data_set(cell_counts_wide,
                                                 cell_metadata=cds_covariates_df,
                                                 gene_metadata=cell_metadata_subset),
                     cds=ccs_cds,
                     cds_coldata=cds_coldata,
                     cds_reduced_dims=cds_reducedDims,
                     info=SimpleList(sample_group=sample_group,
                                     cell_group=cell_group,
                                     norm_method = norm_method))
  
  #
  # PLNmodels::prepare_data returns (1) a matrix of cell abundances,
  # which were calculate in new_cell_count_set() where rows are
  # sample groups and the columns are cell groups, (2) covariates,
  # where is a copy of colData(cds), and (3) offsets, which are
  # calculated by PLNmodels::prepare_data.
  if (norm_method == "size_factors") {
    if (!is.null(size_factors)) {
      assertthat::assert_that(
        tryCatch(expr = identical(sort(colnames(ccs)), sort(names(size_factors))),
                 error = function(e) FALSE),
        msg = "Argument size factor names must match ccs column names.")
      
      pln_data <- PLNmodels::prepare_data(counts = counts(ccs) + pseudocount,
                                          covariates = colData(ccs) %>% as.data.frame,
                                          offset = size_factors)
    } else {
      pln_data <- PLNmodels::prepare_data(counts = counts(ccs) + pseudocount,
                                          covariates = colData(ccs) %>% as.data.frame,
                                          offset = monocle3::size_factors(ccs))
    }
  } else if (norm_method == "RLE") {
    pln_data <- PLNmodels::prepare_data(counts = counts(ccs),
                                        covariates = colData(ccs) %>% as.data.frame,
                                        offset = norm_method,
                                        type="poscounts")
  } else {
    pln_data <- PLNmodels::prepare_data(counts = counts(ccs) + pseudocount,
                                        covariates = colData(ccs) %>% as.data.frame,
                                        offset = norm_method)
    
    if (norm_method == "none") {
      pln_data$Offset = 1
    }
  }
  
  if (norm_method != "size_factors")
    colData(ccs)$Size_Factor = pln_data$Offset
  
  if (keep_cds == FALSE)
    ccs@cds = new_cell_data_set(empty_sparse_matrix(format="C"))
  
  
  # if (!is.null(cell_metadata)) {
  #   assertthat::assert_that(!is.null(row.names(cell_metadata)) &
  #                             all(row.names(cell_metadata) == colnames(expression_data)),
  #                           msg = paste("row.names of cell_metadata must be equal to colnames of",
  #                                       "expression_data"))
  # }
  #
  # if (!is.null(gene_metadata)) {
  #   assertthat::assert_that(!is.null(row.names(gene_metadata)) & all(
  #     row.names(gene_metadata) == row.names(expression_data)),
  #     msg = paste("row.names of gene_metadata must be equal to row.names of",
  #                 "expression_data"))
  # }
  #
  # if (is.null(cell_metadata)) {
  #   cell_metadata <- data.frame(cell = colnames(expression_data),
  #                               row.names = colnames(expression_data))
  # }
  #
  
  # Notes:
  #   o  ccs_cds has the original column names whereas coldata_df has
  #      several renamed columns.
  #   o  coldata_df has all rows
  #   o  ccs_cds has rows filtered by thresholds
  ccs@metadata[["cell_group_assignments"]] = coldata_df %>% dplyr::select(group_id, sample, cell_group) %>% as.data.frame
  ccs@metadata[["cell_group_assignments"]] = ccs@metadata[["cell_group_assignments"]] %>% filter(!cell_group %in% removed_cell_states)
  row.names(ccs@metadata[["cell_group_assignments"]]) = colnames(ccs_cds)
  
  return (ccs)
}




#' Compute a pseudobulk expression matrix for a ccs
#' @export
#' @noRd
pseudobulk_ccs_for_states <- function(ccs, state_col=NULL, collapse_samples=FALSE){
  
  if (is.null(state_col)){
    cell_group_df = tibble::rownames_to_column(ccs@metadata[["cell_group_assignments"]])
    if (collapse_samples)
      cell_group_df = cell_group_df %>% mutate(group_id = cell_group)
    cell_group_df = cell_group_df %>%
      dplyr::mutate(pseudobulk_id = paste(group_id, "cell_group", sep="_")) %>% dplyr::select(rowname, pseudobulk_id, cell_group)
    agg_coldata = cell_group_df %>%
      dplyr::group_by(pseudobulk_id, cell_group) %>%
      dplyr::summarize(num_cells_in_group = n()) %>%
      as.data.frame
    #%>% select(rowname, cell_group)
  }else{
    cell_group_df = tibble::rownames_to_column(ccs@metadata[["cell_group_assignments"]])
    cds_group_df = colData(ccs@cds) %>%
      as.data.frame %>% tibble::rownames_to_column() %>% dplyr::select(rowname, !!sym(state_col))
    cell_group_df = left_join(cell_group_df, cds_group_df, by=c("rowname"))
    if (collapse_samples)
      cell_group_df = cell_group_df %>% mutate(group_id = !!sym(state_col))
    cell_group_df = cell_group_df %>%
      dplyr::mutate(pseudobulk_id = paste(group_id, !!sym(state_col), sep="_")) %>% dplyr::select(rowname, pseudobulk_id, !!sym(state_col))
    agg_coldata = cell_group_df %>%
      dplyr::group_by(pseudobulk_id, !!sym(state_col)) %>%
      dplyr::summarize(num_cells_in_group = n()) %>%
      as.data.frame
  }
  
  agg_expr_mat = monocle3::aggregate_gene_expression(ccs@cds,
                                                     cell_group_df=cell_group_df,
                                                     norm_method="size_only",
                                                     scale_agg_values = FALSE,
                                                     pseudocount=0,
                                                     cell_agg_fun="mean")
  
  agg_expr_mat = agg_expr_mat[,agg_coldata$pseudobulk_id]
  
  row.names(agg_coldata) = agg_coldata$pseudobulk_id
  agg_coldata = agg_coldata[colnames(agg_expr_mat),]
  
  pseudobulk_cds = new_cell_data_set(agg_expr_mat, cell_metadata = agg_coldata, rowData(ccs@cds) %>% as.data.frame)
  pseudobulk_cds = estimate_size_factors(pseudobulk_cds, round_exprs = FALSE)
  return(pseudobulk_cds)
}

#' add metadata to pb_cds from cds
#' @export
add_covariate <- function(ccs, pb_cds, covariate) {
  
  assertthat::assert_that(
    tryCatch(expr = covariate %in% colnames(colData(ccs@cds)),
             error = function(e) FALSE),
    msg = paste0(covariate, " not in colnames"))
  
  assertthat::assert_that(
    tryCatch(expr = !covariate %in% colnames(colData(pb_cds)),
             error = function(e) FALSE),
    msg = paste0(covariate," already in colnames"))
  
  group_to_covariate = colData(ccs@cds) %>%
    as.data.frame %>%
    select(group_id, all_of(covariate)) %>%
    distinct()
  
  pb_coldata = colData(pb_cds) %>%
    as.data.frame %>%
    mutate(group_id = gsub("_cell_group", "", pseudobulk_id)) %>%
    left_join(group_to_covariate, by = "group_id")
  
  colData(pb_cds)[[covariate]] =  pb_coldata[[covariate]]
  
  return(pb_cds)
}



# fit_models(pb_cds,
#            model_formula_str, 
#            weights=colData(pb_cds)$num_cells_in_group)