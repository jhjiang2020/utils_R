#' Fast Wilcoxon rank sum test and auROC 
#' 
#' Adapted from https://github.com/immunogenomics/presto
#' 
#' Computes auROC and Wilcoxon p-value based on Gaussian approximation. 
#' Inputs can be 
#' \itemize{
#' \item Dense matrix or data.frame
#' \item Sparse matrix, such as dgCMatrix
#' \item Seurat V3 object
#' \item SingleCellExperiment object
#' }
#' For detailed examples, consult the presto vignette. 
#' 
#' @param X A feature-by-sample matrix, Seurat object, or SingleCellExperiment
#'  object
#' @param y vector of group labels. 
#' @param groups_use (optional) which groups from y vector to test. 
#' @param group_by (Seurat & SCE) name of groups variable ('e.g. Cluster').
#' @param assay (Seurat & SCE) name of feature matrix slot (e.g. 'data' or
#'  'logcounts'). 
#' @param seurat_assay (Seurat) name of Seurat Assay (e.g. 'RNA'). 
#' @param verbose boolean, TRUE for warnings and messages. 
#' @param ... input specific parameters. 
#' 
#' @examples
#' 
#' data(exprs)
#' data(y)
#' 
#' ## on a dense matrix
#' head(wilcoxauc(exprs, y))
#' 
#' ## with only some groups
#' head(wilcoxauc(exprs, y, c('A', 'B')))
#' 
#' ## on a sparse matrix
#' exprs_sparse <- as(exprs, 'dgCMatrix')
#' head(wilcoxauc(exprs_sparse, y))
#' 
#' ## on a Seurat V3 object
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'     pkg_version <- packageVersion('Seurat')
#'     if (pkg_version >= "3.0" & pkg_version < "4.0") {
#'         data(object_seurat)
#'         head(wilcoxauc(object_seurat, 'cell_type'))
#'     }
#' }
#' 
#' ## on a SingleCellExperiment object
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'     data(object_sce)
#'     head(wilcoxauc(object_sce, 'cell_type'))
#' }
#' 
#' @return table with the following columns: 
#' \itemize{
#' \item \strong{feature} - feature name (e.g. gene name).
#' \item \strong{group} - group name.
#' \item \strong{avgExpr} - mean value of feature in group. 
#' \item \strong{logFC} - log fold change between observations in group vs out.
#' \item \strong{statistic} - Wilcoxon rank sum U statistic. 
#' \item \strong{auc} - area under the receiver operator curve. 
#' \item \strong{pval} - nominal p value. 
#' \item \strong{padj} - Benjamini-Hochberg adjusted p value. 
#' \item \strong{pct_in} - Percent of observations in the group with non-zero
#' feature value. 
#' \item \strong{pct_out} - Percent of observations out of the group with 
#' non-zero feature value. 
#' }
#' @export 

require(Matrix)
wilcoxauc <- function(X, ...) {
  UseMethod('wilcoxauc')
}

#' @rdname wilcoxauc
#' @export
wilcoxauc.seurat <- function(X, ...) {
  stop('wilcoxauc only implemented for Seurat Version 3, please upgrade to 
        run.')
}

#' @rdname wilcoxauc
#' @export
wilcoxauc.Seurat <- function(
    X, group_by=NULL, assay='data', groups_use=NULL, seurat_assay='RNA', ...
) {
  X_matrix <- Seurat::GetAssayData(X, assay=seurat_assay, slot=assay)
  if (is.null(group_by)) {
    y <- Seurat::Idents(X)
  } else {
    y <- Seurat::FetchData(X, group_by) %>% unlist %>% as.character()        
  }
  wilcoxauc(X_matrix, y, groups_use)
}

#' @rdname wilcoxauc
#' @export
wilcoxauc.SingleCellExperiment <- function(
    X, group_by=NULL, assay=NULL, groups_use=NULL, ...
) {
  if (is.null(group_by)) {
    stop('Must specify group_by with SingleCellExperiment')
  } else if (!group_by %in% names(SummarizedExperiment::colData(X))) {
    stop('group_by value is not defined in colData.')
  }
  y <- SummarizedExperiment::colData(X)[[group_by]]
  
  if (is.null(assay)) {
    logcounts <- SingleCellExperiment::logcounts
    standard_assays <- c(
      'normcounts', 'logcounts', 'cpm', 'tpm',
      'weights', 'counts')
    standard_assays <- factor(standard_assays, standard_assays)
    available_assays <- names(SummarizedExperiment::assays(X))
    available_assays <- intersect(standard_assays, available_assays)
    if (length(available_assays) == 0) {
      stop('No assays in SingleCellExperiment object')
    } else {
      assay <- available_assays[1]
    }
  }
  
  X_matrix <- eval(call(assay, X))
  wilcoxauc(X_matrix, y, groups_use)
}

#' @rdname wilcoxauc
#' @export
wilcoxauc.default <- function(X, y, groups_use=NULL, verbose=TRUE, ...) {
  ## Check and possibly correct input values
  if (is(X, 'dgeMatrix')) X <- as.matrix(X)
  if (is(X, 'data.frame')) X <- as.matrix(X)
  # if (is(X, 'DataFrame')) X <- as.matrix(X)
  # if (is(X, 'data.table')) X <- as.matrix(X)
  if (is(X, 'dgTMatrix')) X <- as(X, 'dgCMatrix')
  if (is(X, 'TsparseMatrix')) X <- as(X, 'dgCMatrix')
  if (ncol(X) != length(y)) stop("number of columns of X does not
                                match length of y")
  if (!is.null(groups_use)) {
    idx_use <- which(y %in% intersect(groups_use, y))
    y <- y[idx_use]
    X <- X[, idx_use]
  }
  
  
  y <- factor(y)
  idx_use <- which(!is.na(y))
  if (length(idx_use) < length(y)) {
    y <- y[idx_use]
    X <- X[, idx_use]
    if (verbose) 
      message('Removing NA values from labels')        
  }
  
  group.size <- as.numeric(table(y))
  if (length(group.size[group.size > 0]) < 2) {
    stop('Must have at least 2 groups defined.')
  }
  
  #     features_use <- which(apply(!is.na(X), 1, all))
  #     if (verbose & length(features_use) < nrow(X)) {
  #         message('Removing features with NA values')
  #     }
  #     X <- X[features_use, ]
  if (is.null(row.names(X))) {
    row.names(X) <- paste0('Feature', seq_len(nrow(X)))
  }
  
  ## Compute primary statistics
  group.size <- as.numeric(table(y))
  n1n2 <- group.size * (ncol(X) - group.size)

  if (is(X, 'dgCMatrix')) {
    rank_res <- rank_matrix(Matrix::t(X))    
    cohen_d <- compute_cohend(Matrix::t(X), y, group.size)
  } else {
    rank_res <- rank_matrix(X)
    cohen_d <- compute_cohend(X, y, group.size)
  }
  
  ustat <- compute_ustat(rank_res$X_ranked, y, n1n2, group.size) 
  auc <- t(ustat / n1n2)
  pvals <- compute_pval(ustat, rank_res$ties, ncol(X), n1n2) 
  fdr <- apply(pvals, 2, function(x) p.adjust(x, 'BH'))
  

  
  ### Auxiliary Statistics (AvgExpr, PctIn, LFC, etc)
  group_sums <- sumGroups(X, y, 1)
  group_nnz <- nnzeroGroups(X, y, 1)
  group_pct <- sweep(group_nnz, 1, as.numeric(table(y)), "/") %>% t()
  group_pct_out <- -group_nnz %>% 
    sweep(2, colSums(group_nnz) , "+") %>% 
    sweep(1, as.numeric(length(y) - table(y)), "/") %>% t()
  group_means <- sweep(group_sums, 1, as.numeric(table(y)), "/") %>% t()
  cs <- colSums(group_sums)
  gs <- as.numeric(table(y))
  lfc <- Reduce(cbind, lapply(seq_len(length(levels(y))), function(g) {
    group_means[, g] - ((cs - group_sums[g, ]) / (length(y) - gs[g]))
  }))
  
  res_list <- list(auc = auc, 
                   pval = pvals,
                   padj = fdr, 
                   pct_in = 100 * group_pct, 
                   pct_out = 100 * group_pct_out,
                   avgExpr = group_means, 
                   statistic = t(ustat),
                   logFC = lfc,
                   cohen_d = t(cohen_d))
  return(tidy_results(res_list, row.names(X), levels(y)))
}


#' Get top n markers from wilcoxauc
#' 
#' Useful summary of the most distinguishing features in each group. 
#' 
#' @param res table returned by wilcoxauc() function. 
#' @param n number of markers to find for each. 
#' @param auc_min filter features with auc < auc_min. 
#' @param pval_max filter features with pval > pval_max.
#' @param padj_max  filter features with padj > padj_max.
#' @param pct_in_min Minimum percent (0-100) of observations with non-zero 
#' entries in group.
#' @param pct_out_max Maximum percent (0-100) of observations with non-zero 
#' entries out of group.
#' @examples
#' 
#' data(exprs)
#' data(y)
#' 
#' ## first, run wilcoxauc
#' res <- wilcoxauc(exprs, y)
#' 
#' ## top 10 markers for each group
#' ## filter for nominally significant (p<0.05) and over-expressed (auc>0.5)
#' top_markers(res, 10, auc_min = 0.5, pval_max = 0.05)
#' 
#' @return table with the top n markers for each cluster. 
#' @export 
top_markers <- function(res, n=10, auc_min=0, pval_max=1, padj_max=1,
                        pct_in_min=0, pct_out_max=100) {
  res %>% 
    dplyr::filter(
      pval <= pval_max & 
        padj <= padj_max &
        auc >= auc_min & 
        pct_in >= pct_in_min &
        pct_out <= pct_out_max
    ) %>%
    dplyr::group_by(group) %>%
    dplyr::top_n(n = n, wt = auc) %>% 
    dplyr::mutate(rank = rank(-auc, ties.method = 'random')) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(feature, group, rank) %>% 
    tidyr::spread(group, feature, fill = NA)
}


#' rank_matrix
#' 
#' Utility function to rank columns of matrix
#' 
#' @param X feature by observation matrix. 
#' 
#' @examples
#' 
#' data(exprs)
#' rank_res <- rank_matrix(exprs)
#' 
#' @return List with 2 items
#' \itemize{
#' \item X_ranked - matrix of entry ranks
#' \item ties - list of tied group sizes
#' }
#' @export 
rank_matrix <- function(X) {
  UseMethod('rank_matrix')
}

#' @rdname rank_matrix
#' @export
rank_matrix.dgCMatrix <- function(X) {
  Xr <- Matrix(X, sparse = TRUE)
  ties <- cpp_rank_matrix_dgc(Xr@x, Xr@p, nrow(Xr), ncol(Xr))
  return(list(X_ranked = Xr, ties = ties))
}

#' @rdname rank_matrix
#' @export
rank_matrix.matrix <- function(X) {
  cpp_rank_matrix_dense(X)
}

compute_ustat <- function(Xr, cols, n1n2, group.size) {
  grs <- sumGroups(Xr, cols)
  if (is(Xr, 'dgCMatrix')) {
    gnz <- (group.size - nnzeroGroups(Xr, cols))
    zero.ranks <- (nrow(Xr) - diff(Xr@p) + 1) / 2
    ustat <- t((t(gnz) * zero.ranks)) + grs - group.size *
      (group.size + 1 ) / 2        
  } else {
    ustat <- grs - group.size * (group.size + 1 ) / 2
  }
  return(ustat)
}

# cohen's d = (m1-m2)/sd_pooled (assumed two normal distributions)
compute_cohend <- function(Xr, cols, group.size) {
  m1 <- sumGroups(Xr, cols)/group.size
  m2 <- -sweep(sumGroups(Xr, cols), MARGIN = 2, STATS = colSums(sumGroups(Xr, cols)))/ (sum(group.size)-group.size)
  cohen_d <- (m1-m2)/sqrt(MatrixGenerics::colVars(Xr))
  return(cohen_d)
}

compute_pval <- function(ustat, ties, N, n1n2) {
  z <- ustat - .5 * n1n2
  z <- z - sign(z) * .5
  .x1 <- N ^ 3 - N
  .x2 <- 1 / (12 * (N^2 - N))
  rhs <- lapply(ties, function(tvals) {
    (.x1 - sum(tvals ^ 3 - tvals)) * .x2
  }) %>% unlist
  usigma <- sqrt(matrix(n1n2, ncol = 1) %*% matrix(rhs, nrow = 1))
  z <- t(z / usigma)
  
  pvals <- matrix(2 * pnorm(-abs(as.numeric(z))), ncol = ncol(z))
  return(pvals)
}

#' sumGroups
#' 
#' Utility function to sum over group labels
#' 
#' @param X matrix
#' @param y group labels
#' @param MARGIN whether observations are rows (=2) or columns (=1)
#' 
#' @examples
#' 
#' data(exprs)
#' data(y)
#' sumGroups_res <- sumGroups(exprs, y, 1)
#' sumGroups_res <- sumGroups(t(exprs), y, 2)
#' 
#' @return Matrix of groups by features
#' @export 
sumGroups <- function(X, y, MARGIN=2) {
  if (MARGIN == 2 & nrow(X) != length(y)) {
    stop('wrong dims')
  } else if (MARGIN == 1 & ncol(X) != length(y)) {
    stop('wrong dims') 
  }
  UseMethod('sumGroups')
}

#' @rdname sumGroups
#' @export
sumGroups.dgCMatrix <- function(X, y, MARGIN=2) {
  if (MARGIN == 1) {
    cpp_sumGroups_dgc_T(X@x, X@p, X@i, ncol(X), nrow(X), as.integer(y) - 1,
                        length(unique(y)))        
  } else {
    cpp_sumGroups_dgc(X@x, X@p, X@i, ncol(X), as.integer(y) - 1,
                      length(unique(y)))
  }
}

#' @rdname sumGroups
#' @export
sumGroups.matrix <- function(X, y, MARGIN=2) {
  if (MARGIN == 1) {
    cpp_sumGroups_dense_T(X, as.integer(y) - 1, length(unique(y)))        
  } else {
    cpp_sumGroups_dense(X, as.integer(y) - 1, length(unique(y)))
  }
}



tidy_results <- function(wide_res, features, groups) {
  res <- Reduce(cbind, lapply(wide_res, as.numeric)) %>% data.frame() 
  colnames(res) <- names(wide_res)
  res$feature <- rep(features, times = length(groups))
  res$group <- rep(groups, each = length(features))
  res %>% dplyr::select(
    feature, 
    group, 
    avgExpr, 
    logFC, 
    cohen_d,
    statistic, 
    auc, 
    pval, 
    padj, 
    pct_in, 
    pct_out
  )
}


#' nnzeroGroups
#' 
#' Utility function to compute number of zeros-per-feature within group
#' 
#' @param X matrix
#' @param y group labels
#' @param MARGIN whether observations are rows (=2) or columns (=1)
#' 
#' @examples
#' 
#' data(exprs)
#' data(y)
#' nnz_res <- nnzeroGroups(exprs, y, 1)
#' nnz_res <- nnzeroGroups(t(exprs), y, 2)
#' 
#' @return Matrix of groups by features
#' @export 
nnzeroGroups <- function(X, y, MARGIN=2) {
  if (MARGIN == 2 & nrow(X) != length(y)) {
    stop('wrong dims')
  } else if (MARGIN == 1 & ncol(X) != length(y)) {
    stop('wrong dims')        
  }
  UseMethod('nnzeroGroups')
}

#' @rdname nnzeroGroups
#' @export
nnzeroGroups.dgCMatrix <- function(X, y, MARGIN=2) {
  if (MARGIN == 1) {
    cpp_nnzeroGroups_dgc_T(X@p, X@i, ncol(X), nrow(X), as.integer(y) - 1,
                           length(unique(y)))        
  } else {
    cpp_nnzeroGroups_dgc(X@p, X@i, ncol(X), as.integer(y) - 1,
                         length(unique(y)))
  }
}

#' @rdname nnzeroGroups
#' @export
nnzeroGroups.matrix <- function(X, y, MARGIN=2) {
  if (MARGIN == 1) {        
    cpp_nnzeroGroups_dense_T(X, as.integer(y) - 1, length(unique(y)))        
  } else {
    cpp_nnzeroGroups_dense(X, as.integer(y) - 1, length(unique(y)))
  }
}
