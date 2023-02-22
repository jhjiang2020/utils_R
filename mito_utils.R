# customised function based on Signac::AlleleFreq
#
AlleleSetFreq <- function(object, variants, return_DP=FALSE, ...) {
  variants <- unique(x = variants)
  # Access meta data for the counts
  meta_row_mat <- as.data.frame(
    x = stringi::stri_split_fixed(
      str = rownames(x = object),
      pattern = "-",
      simplify = TRUE
    ), stringsAsFactors = TRUE
  )
  colnames(meta_row_mat) = c("letter", "position", "strand")
  
  # Access meta data for the variants
  variant_df <- data.frame(
    variant = variants,
    position = factor(
      x = substr(
        x = variants,
        start = 1,
        stop = nchar(x = variants) - 3),
      levels = levels(x = meta_row_mat$position)
    ),
    ref = factor(
      x = substr(
        x = variants,
        start = nchar(x = variants) - 2,
        stop = nchar(x = variants) - 2),
      levels = levels(x = meta_row_mat$letter)
    ),
    alt = factor(
      x = substr(
        x = variants,
        start = nchar(x = variants),
        stop = nchar(x = variants)),
      levels = levels(x = meta_row_mat$letter)
    )
  )
  
  # Numerator counts
  # Get the forward and reverse strands for the matching the position / letter
  # for the alternate allele
  ref_letter <-  paste0(meta_row_mat$position, meta_row_mat$letter)
  alt_letter <- paste0(variant_df$position, variant_df$alt)
  idx_numerator <- lapply(
    X = alt_letter, FUN = function(x) {
      which(ref_letter == x)
    }
  )
  fwd_half_idx <- sapply(X = idx_numerator, FUN = `[[`, 1)
  rev_half_idx <- sapply(X = idx_numerator, FUN = `[[`, 2)
  
  # verify that the object is behaving like we expect
  if (!all.equal(
    target = meta_row_mat[fwd_half_idx, 2],
    current = meta_row_mat[rev_half_idx, 2]
  )) {
    stop("Variant count matrix does not have the required structure")
  }
  numerator_counts <- object[fwd_half_idx, ] + object[rev_half_idx, ]
  rownames(x = numerator_counts) <- variants
  
  # Same idea for the denominator but use all counts at each position
  denom_counts <- sapply(X = variant_df$position, FUN = function(x) {
    idx <- which(meta_row_mat$position == x)
    total_coverage <- Matrix::colSums(x = object[idx, ])
    return(total_coverage)
  })
  denom_counts <- t(x = denom_counts)
  rownames(x = denom_counts) <- variant_df$variant
  
  # final frequency as the sum of total numerator / the sum of total denominator
  alleleset_freq <- Matrix::colSums(numerator_counts) / Matrix::colSums(denom_counts)
  
  # Set NaN value due to 0 total counts to 0
  alleleset_freq[is.nan(x = alleleset_freq)] <- 0
  
  names(alleleset_freq) <- colnames(object)
  
  # return total sequencing depth aggregating all variants
  if(return_DP){
    return(data.frame("dp" = Matrix::colSums(denom_counts),
                      "freq" = alleleset_freq
                      ))
  }
  return(alleleset_freq)
}

# Run correlation test in parallel 
runCor <- function(rnaMat, obs, cells_select=NULL, nCores=1, chunksize=1, method="spearman"){
  cat("Testing ",nrow(rnaMat)," Genes\n")
  stopifnot(all(colnames(rnaMat)==names(obs)))
  names(obs) <- colnames(rnaMat)
  if(!is.null(cells_select)){
    cat("Testing in ",length(cells_select)," cells\n")
    stopifnot(all(cells_select %in% colnames(rnaMat)))
    rnaMat <- rnaMat[,cells_select]
    obs <- obs[cells_select]
  }
  library(doParallel)
  if(nCores > 1)
    message("Running correlation using ",nCores," cores ..\n")
  opts <- list()
  pb <- txtProgressBar(min = 0, max = nrow(rnaMat), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  time_elapsed <- Sys.time()
  cl <- parallel::makeCluster(nCores)
  doSNOW::registerDoSNOW(cl)
  
  
  mCortest.list <- foreach(g=rownames(rnaMat),
                           .options.snow = opts, 
                           #"BuenRTools"
                           .packages = c( "dplyr","Matrix")) %dopar%   {
                             # Correlate smoothed Gene expression with HTP percentage, with spearman
                             corr.r <- cor(rnaMat[g,], obs, method = method)
                             rownames(corr.r) <- g
                             return(corr.r)
                           }
  cat("Finished!\n")
  cat("Merging results ..\n")
  # Merge and save table for downstream filtering and plotting (network)
  GeneCor.d <- do.call('rbind',mCortest.list)
  return(GeneCor.d)
}
