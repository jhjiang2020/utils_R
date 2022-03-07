require(dplyr)
format_GWAS_SNP <- function(catalogfile, p_value_threshold){
  chrom <- chromEnd <- SNP <- marker <- P <- trait <- NULL
  gwas.hg38 <- readr::read_tsv(catalogfile, col_names = FALSE, skip = 1)
  if (ncol(gwas.hg38) != 23) stop("The wrong number of columns are present in the
                                file - are you sure this is the correct file?")
  
  colnames(gwas.hg38) <- c("bin", "chrom", "chromStart", "chromEnd", "SNP",
                           "pubMedID", "author", "pubDate", "journal", "title",
                           "trait", "initSample", "replSample", "region",
                           "genes", "riskAllele", "riskAlFreq", "P",
                           "pValueDesc", "orOrBeta", "ci95", "platform",
                           "cnv")
  
  gwas.hg38.formatted <- gwas.hg38 %>%
    dplyr::mutate(marker = paste(gsub("chr", "", chrom), chromEnd, sep = ":")) %>%   ## 
    dplyr::select(SNP, marker, P, trait) %>%
    dplyr::arrange(P) %>%
    dplyr::filter(P < p_value_threshold)
  
  return(gwas.hg38.formatted)
}

filter_GWAS_SNP <- function(formatted.gwas, filter){
  # filter format as "termA|termB" 
  filter.criteria <- gsub (" \\| ", "\\|", filter)
  filter.criteria <- trimws(filter.criteria)
  return(formatted.gwas[grepl(filter.criteria,formatted.gwas$trait,ignore.case = TRUE), ])
}

prune_GWAS_SNP <- function(plink = NULL, snps, genotypeData) {
  
  # Check that plink command works
  tryCatch({
    null <- system(command = paste(plink, "--help"), intern = TRUE)
  }, error = function(cond) {
    message(paste("Software does not seem to exist:", plink))
    message("Here's the original error message:")
    message(cond)
  } )
  
  if (is.data.frame(snps)) {
    tmp_name <- tempfile(pattern = "tmp", tmpdir = tempdir(), fileext = ".txt")
    readr::write_tsv(snps, file = tmp_name, col_names = TRUE) # write_tsv is faster than the base R write.table function
  }else if(is.character(snps)){
    warning("No P value provided, disabling filtering......\n")
    tmp_name <- tempfile(pattern = "tmp", tmpdir = tempdir(), fileext = ".txt")
    readr::write_tsv(data.frame("SNP" = snps, "P" = 5e-8), file = tmp_name, col_names = T)
  }else{
    stop("Input SNP must be a dataframe or a string!")
  }
  
  #gtdf <- file.path(genotypeData, 
  #"ALL.chr_merged.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes")
  gtdf <- genotypeData
  out_file <- tempfile(pattern = "clumpedSNPs", tmpdir = tempdir())
  code.clumping <- sprintf(
    "%s --bfile %s --clump-p1 5e-8 --clump-kb 500 --clump-r2 0.1 --clump %s --out %s",
    plink, gtdf, tmp_name, out_file)
  system(code.clumping)
  snps.clumped <- data.table::fread(paste0(out_file,".clumped"),
                                    header = T, stringsAsFactors = F)
  unlink(tmp_name)
  unlink(list.files(pattern = paste0(out_file, "*")))
  return(snps.clumped)
}

getld_GWAS_SNP <- function(plink, genotypeData, snps, r2 = 0.8, return_clump = FALSE) {
  
  CHR_A <- BP_A <- CHR_B <- BP_B <- R2 <- SNP_B <- comb <- indexSNP <- ldPos <- NULL
  
  if (is.data.frame(snps)) {
    tmp_name <- tempfile(pattern = "tmp", tmpdir = tempdir(), fileext = ".txt")
    readr::write_tsv(snps, file = tmp_name, col_names = TRUE) # write_tsv is faster than the base R write.table function
  }else if(is.character(snps)){
      if(any(duplicated(snps))){
        warning("Removing duplicated rsids......\n")
        snps <- unique(snps)
      }
    tmp_name <- tempfile(pattern = "tmp", tmpdir = tempdir(), fileext = ".txt")
    readr::write_tsv(data.frame("SNP" = snps), file = tmp_name, col_names = T)
  }else{
    stop("Input SNP must be a dataframe or a string!")
  }
  
  out_prefix <- paste0("ldpartners_", gsub("\\.", "", r2))
  out <- tempfile(pattern = out_prefix, tmpdir = tempdir() )
  
  gtdf <- genotypeData
  
  #file.path(genotypeData, 
  #                "ALL.chr_merged.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes") 
  
  ## call plink to get the ld of input snps
  code.ld <- sprintf(
    "%s --bfile %s --r2 --ld-snp-list %s --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 %s --out %s",
    plink, gtdf, tmp_name, r2, out)
  system(code.ld)
  ld.partners <- data.table::fread(paste0(out, ".ld"), header = T,
                                   stringsAsFactors = FALSE)
  ld.partners.r2 <- ld.partners
  ld.partners.r2[, c("indexPos", "ldPos") := {
    list(paste(CHR_A, BP_A, sep = ":"),
         paste(CHR_B, BP_B, sep = ":"))
    }
  ]
  ld.partners.r2 <- ld.partners.r2[, c("SNP_A", "indexPos", "SNP_B", "ldPos", "R2")]
  oldnames <- colnames(ld.partners.r2)
  newnames <- c("indexSNP", "indexPos", "ldSNP", "ldPos", "r2")
  data.table::setnames(ld.partners.r2, old = oldnames, new = newnames)
  
  if(return_clump){ 
    warning("Make sure the SNP list is pruned before running this command......\n")
    clump <- snps$SP2
    clump <- lapply(clump, function(x) {
        x <- unlist(strsplit(x, ","))
        x <- gsub("\\(1\\)", "", x)
      }
    )
    ld <- ld.partners[, 6:7]
    clump.ld <- lapply(clump, function(x) ld[ld$SNP_B %in% x, ])
    clump.ld <- lapply(clump.ld, function(x) {
      x %>% dplyr::arrange(desc(R2)) %>%
        dplyr::mutate(comb = paste0(SNP_B, " (", round(R2, 2), ")")) %>%
        dplyr::select(comb)
    })
    clump.ld <- lapply(clump.ld, function(x) {
      x <- do.call(c, x)
      x <- paste(x, collapse = ", ")
    })
    snps$SP2_ld <- unlist(clump.ld)
    results <- ld.partners.r2 %>% dplyr::group_by(indexSNP) %>%
      dplyr::summarise(locus_up = min(ldPos), locus_down = max(ldPos))
    
    results <- merge(independent_snps, results, by.x = "SNP", by.y = "indexSNP",
                     keep.all = T)
    results.snps <- results[ , c(1, 2, 4, 5, 13, 14, 15)]
    oldnames <- colnames(results.snps)
    newnames <- c('snp_name', 'chr', 'pos','pvalue',
                  'plink_ld_partners', 'locus_upstream_boundary',
                  'locus_downstream_boundary')
    data.table::setnames(results.snps, oldnames, newnames)
    
    unlink(c(tmp_name, list.files(pattern = paste0(out, "*"))))
    return(list(results.snps, ld.partners.r2))
  }
  
  unlink(c(tmp_name, list.files(pattern = paste0(out, "*"))))
  return(ld.partners.r2)
}

format_goshifter_SNPs <- function(snpfinallist){
  ldpos <- snpfinallist[[2]]$"ldPos" %>% 
    strsplit(split = ":") %>% 
    unlist %>% 
    matrix(nrow=dim(snpfinallist[[2]])[1],ncol=2,byrow = T)
  snp_map <- data.frame("SNP"=snpfinallist[[2]]$"ldSNP", "Chrom"=paste("chr",ldpos[,1],sep = ""), "BP"=ldpos[,2])
  return(snp_map)
}

liftover_GWAS_SNP <- function(SNPs, biomaRt_matrix = NULL, local = FALSE){
  if(local){
    return(liftover_GWAS_SNP_local(SNPs))
  }
  if(is.null(biomaRt_matrix)){
    require(biomaRt)
    ensembl <- useEnsembl("snp",dataset = "hsapiens_snp")
    if(is.data.frame(SNPs)){
      stopifnot("ldSNP" %in% colnames(SNPs))
      rsids <- SNPs$ldSNP
    }else if(is.character(SNPs)){
      rsids <- SNPs
    }else{
      stop("Input SNP must be a dataframe or a string vector!")
    }
    #get genomic position
    SNPs <- getBM(attributes=c("refsnp_id",
                               "chr_name",
                               "chrom_start",
                               "chrom_end"),
                  filters ="snp_filter", 
                  values =rsids, 
                  mart = ensembl, uniqueRows=TRUE)
  }else{
    SNPs <- filter(biomaRt_matrix, refsnp_id %in% SNPs)
  }
  SNPs <- SNPs[SNPs$chr_name %in% 1:22,]
  SNPs$chr_name <- paste("chr", SNPs$chr_name, sep = "")
  
  ## This is to make sure the start position always <= end position of any given snps
  select <- SNPs$chrom_start>=SNPs$chrom_end
  SNPs[select,c('chrom_start','chrom_end')] <- SNPs[select,c('chrom_end','chrom_start')]
  
  return(SNPs)
}

liftover_GWAS_SNP_local <- function(SNPs){
  require(SNPlocs.Hsapiens.dbSNP151.GRCh38)
  require(dplyr)
  
  if(is.data.frame(SNPs)){
    stopifnot("ldSNP" %in% colnames(SNPs))
    rsids <- SNPs$ldSNP
  }else if(is.character(SNPs)){
    rsids <- SNPs
  }else{
    stop("Input SNP must be a dataframe or a string vector!")
  }
  
  keep <- gsub(pattern = "(rs)+[0-9]+", replacement = "\\1", rsids) == "rs"
  n_remove = sum(!keep)
  warning(paste0("Removing ", n_remove, " SNPs that do not have rsid...\n" ))
  rsids <- rsids[keep]

  
  SNPs <- snpsById(SNPlocs.Hsapiens.dbSNP151.GRCh38,rsids, ifnotfound="drop") %>% 
    as.data.frame() %>% 
    dplyr::filter(seqnames %in% 1:22) %>% 
    dplyr::mutate(chr_name = paste("chr", seqnames, sep = "")) %>% 
    dplyr::mutate(chrom_start = pos,chrom_end = pos, refsnp_id = RefSNP_id) %>% 
    dplyr::select(refsnp_id, chr_name, chrom_start, chrom_end)
  return(SNPs)
}

snp2gr <- function(SNPs){
  require(GenomicRanges)
  stopifnot(is.data.frame(SNPs))
  stopifnot(ncol(SNPs) == 4)
  colnames(SNPs) <- c("rsid", "chr", "start", "end")
  SNPs <- makeGRangesFromDataFrame(SNPs,
                                  keep.extra.columns=TRUE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field=c("seqnames", "seqname",
                                                   "chromosome", "chrom",
                                                   "chr", "chromosome_name",
                                                   "seqid"),
                                  start.field="start",
                                  end.field=c("end", "stop"),
                                  strand.field="strand",
                                  starts.in.df.are.0based=FALSE)
  return(SNPs)
}

# read_matched_SNPs <- function(){ ### NOT FINISHED!
#   
#   ### matched SNP set should be obtained before getting ld snps
#   query_similar_SNPs <- function(SNPsnap_collection, nperm = 1000, ld.snps) {
#     
#     collection <- ldSNP <- ldPos <- snp_maf <- gene_count <- friends_ld08 <- NULL
#     HGNC_nearest_gene_snpsnap <- dist_nearest_gene_snpsnap <- gene_count <- NULL
#     rowid <- J <- snpID <- end <- chr <- start <- rsID <- NULL
#     
#     cat("\n\nGenerating matched SNPs for permutation testing.  SNPs will be 
#       matched for MAF, number of genes in SNP locus, and LD 'buddies' at 
#       r2 0.8.\n")
#     cat("Number of permutations =", nperm)
#     cat("\nLoading SNP database - this may take a while.  If this takes longer than 2-3 
#       minutes, you may not have sufficient RAM to proceed.\n")
#     load(snpdatabase)
#     cat("\nDatabase loaded.\n")
#     
#     ld.partners.r2.dt <- data.table::data.table(ld.snps[[2]])
#     ld.partners.r2.dt <- lapply(collection, function(x) {
#       tmp <- merge(ld.partners.r2.dt, x, by.x = "ldPos", 
#                    by.y = "snpID", keep.x = T) %>%
#         dplyr::select(ldSNP, ldPos, snp_maf, gene_count, friends_ld08, 
#                       nearest_gene = HGNC_nearest_gene_snpsnap, 
#                       dist_nearest_gene = dist_nearest_gene_snpsnap)
#     })
#     
#     cat("\nCreating sets of matched SNPs for each permutation - this step may be 
#       time-consuming.\n")
#     
#     matched <- vector("list", 49)
#     
#     ################## TO DO: rewrite this as a foreach loop
#     for (i in seq_along(ld.partners.r2.dt)) {
#       x <- ld.partners.r2.dt[[i]]
#       y <- list()
#       for (j in seq_len(nrow(x))) {
#         gc <- as.numeric(x[j, gene_count])
#         ld <- as.numeric(x[j, friends_ld08])
#         s1 <- collection[[i]][as.numeric(gene_count) < 1.3*gc & 
#                                 as.numeric(gene_count) > 0.7*gc]
#         s2 <- s1[as.numeric(friends_ld08) < 1.3*ld & 
#                    as.numeric(friends_ld08) > 0.7*ld]
#         y[[j]] <- s2[sample(nrow(s2), nperm, replace = T), ]
#       }
#       matched[[i]] <- y
#       names(matched)[i] <- names(collection)[i]
#     }
#     
#     matched2 <- unlist(matched, recursive = F)
#     # This bit where we take the ith row from each dataframe in matched2 and rbind
#     # potentially takes forever!  I have finally come up with a data.table solution
#     # that is about 1000x faster!!  
#     
#     # convert to data.table
#     invisible(lapply(matched2, data.table::setattr, name = "class", 
#                      value = c("data.table", "data.frame")))
#     # make one big table
#     bigdata <- data.table::rbindlist(matched2)
#     # generate an index - the number of dataframes will be the number of 
#     # permutations
#     index <- as.character(seq_len(nperm))
#     bigdata[, `:=`(rowid, index)]
#     # set a key based on the row index
#     data.table::setkey(bigdata, rowid)
#     # split on this
#     matched.list <- lapply(index, function(i, j, x) x[i = J(i)], x = bigdata)
#     # Drop the `row id` column
#     invisible(lapply(matched.list, function(x) data.table::set(x, j = "rowid", value = NULL)))
#     # Convert back to data.frame
#     invisible(lapply(matched.list, data.table::setattr, name = 'class', 
#                      value = c('data.frame')))
#     
#     # Convert the matched SNPs into a GRanges object list
#     matched.list <- 
#       lapply(matched.list, function(x) {
#         a <- x %>% dplyr::mutate(chr = paste0("chr", gsub(":.*$", "", snpID)), 
#                                  end = as.numeric(gsub("^.*:", "", snpID)), 
#                                  start = (end - 1)) %>%
#           dplyr::select(chr, start, end, snp = rsID, 
#                         nearest_gene = HGNC_nearest_gene_snpsnap, 
#                         dist_nearest_gene = dist_nearest_gene_snpsnap)
#         a <- a[stats::complete.cases(x), ]
#         return(a)
#       })
#     
#     matched.list.GR <- lapply(matched.list, function(x) {
#       a <- as.data.frame(x)
#       GenomicRanges::makeGRangesFromDataFrame(a,  keep.extra.columns = TRUE, 
#                                               seqnames.field = "chr", start.field = "start",
#                                               end.field = "end")
#     })
#     save(matched.list.GR, file = paste0("matched_SNPs_GR.rda"))
#     return(matched.list.GR)
#   }}


GWASPeakZtest <- function(peakSet, ### peak_set: vector of peak ids you wish to test for motif enrichment
                          bgPeaks, ### bg_peaks: matrix of background peak selection iterations by chromVAR
                          SNPSet, ### SNP_set: a GRange object of GWAS associated SNPs 
                          n_bgs = ncol(bgPeaks), ### optional: number of background peaksets, by default all bgsets in bgPeak matrix
                          return_bg = F, ### optional: if set TRUE, a list containing background overlaps will be returned
                          weights = NULL ### optional: if provided, number of overlaps will be multiplied by this vector. 
) {
  # if no weights provided, enrichment analysis will be performed using # of overlaps
  if(is.null(weights))
    weights <- rep(1, length(peakSet))
  stopifnot(length(weights) == length(peakSet))
  
  # peaks should be in a character vector
  stopifnot(is.character(peakSet))
  
  # for duplicated peaks, only keeping the peak with highest weight
  if(sum(duplicated(peakSet))){
    warning("Removing duplicated peaks from input peak set ..\n")
    df <- data.frame(peaks = peakSet,
                     weights = weights) %>% 
      arrange(peakSet, desc(weights))
    peakSet <- df$peaks[!duplicated(df$peaks)]
    warning("Only keeping the highest weight for each peak ..\n")
    weights <- df$weights[!duplicated(df$peaks)]
  }
  if(ncol(bgPeaks) < n_bgs)
    stop("Backgroud peak matrix must have sufficient number of peakset..\n")
  
  if(!all(peakSet %in% rownames(bgPeaks)))
    stop("One or more of the provided peak indices are out of the background peak set range ..\n")
  
  o <- findOverlaps(SNPSet, peakSet %>% StringToGRanges)
  nLOverlaps <- weights[subjectHits(o)] %>% sum
  
  peaks.ix <- rownames(bgPeaks) %in% peakSet
  # a for loop to test # of bg interations 
  cat("Start", n_bgs, "permutation tests ..\n")
  nBGOverlaps <- c() 
  for( i in 1:n_bgs){
    bgPeaks.ix <- bgPeaks[peaks.ix,i]
    bgPeaks.set <- rownames(bgPeaks)[bgPeaks.ix] %>% StringToGRanges
    o <- findOverlaps(SNPSet, bgPeaks.set)
    nOverlaps <- weights[subjectHits(o)] %>% sum
    nBGOverlaps <- c(nBGOverlaps, nOverlaps)
  }
  OR <- nLOverlaps/mean(nBGOverlaps)
  cat("Calculate Z-score based on background distribution ..\n")
  z_score <- (nLOverlaps - mean(nBGOverlaps)) / sd(nBGOverlaps)
  pval.perm = sum(nBGOverlaps > nLOverlaps)/length(nBGOverlaps)
  
  d <- data.frame(
    # nCADSNPs = length(o),
    Z = z_score,
    pval.z = pnorm(-abs(z_score)), ## one-tailed test ,
    signed.log10p = -log10( pnorm(-abs(z_score))) * sign(z_score), 
    pval.perm = pval.perm,
    signed.log10pperm = -log10(pval.perm) * sign(pval.perm),
    Obs.overlaps = nLOverlaps,
    OR = OR
  )
  if(!return_bg){return(d)}
  else{
    return(list(d,nBGOverlaps))
  }
}

MarkPeakZtest <- function(peakSet, ### peak_set: vector of peak ids you wish to test for motif enrichment
                          bgPeaks, ### bg_peaks: matrix of background peak selection iterations by chromVAR
                          MarkSet, ### Mark_set: a GRange object of marker
                          n_bgs = ncol(bgPeaks), ### optional: number of background peaksets, by default all bgsets in bgPeak matrix
                          return_bg = F, ### optional: if set TRUE, a list containing background overlaps will be returned
                          weights = NULL ### optional: if provided, number of overlaps will be multiplied by this vector. 
) {
  # if no weights provided, enrichment analysis will be performed using # of overlaps
  if(is.null(weights))
    weights <- rep(1, length(peakSet))
  stopifnot(length(weights) == length(peakSet))
  
  # for duplicated peaks, only keeping the peak with highest weight
  if(sum(duplicated(peakSet))){
    warning("Removing duplicated peaks from input peak set ..\n")
    df <- data.frame(peaks = peakSet,
                     weights = weights) %>% 
      arrange(peakSet, desc(weights))
    peakSet <- df$peaks[!duplicated(df$peaks)]
    warning("Only keeping the highest weight for each peak ..\n")
    weights <- df$weights[!duplicated(df$peaks)]
  }
  if(ncol(bgPeaks) < n_bgs)
    stop("Backgroud peak matrix must have sufficient number of peakset..\n")
  
  if(!all(peakSet %in% rownames(bgPeaks)))
    stop("One or more of the provided peak indices are out of the background peak set range ..\n")
  
  
  o <- findOverlaps(MarkSet, peakSet %>% StringToGRanges)
  nLOverlaps <- weights[subjectHits(o)] %>% sum
  
  peaks.ix <- rownames(bgPeaks) %in% peakSet
  # a for loop to test # of bg interations 
  cat("Start", n_bgs, "permutation tests ..\n")
  nBGOverlaps <- c() 
  for( i in 1:n_bgs){
    bgPeaks.ix <- bgPeaks[peaks.ix,i]
    bgPeaks.set <- rownames(bgPeaks)[bgPeaks.ix] %>% StringToGRanges
    o <- findOverlaps(MarkSet, bgPeaks.set)
    nOverlaps <- weights[subjectHits(o)] %>% sum
    nBGOverlaps <- c(nBGOverlaps, nOverlaps)
  }
  
  
  OR <- nLOverlaps/mean(nBGOverlaps)
  cat("Calculate Z-score based on background distribution ..\n")
  z_score <- (nLOverlaps - mean(nBGOverlaps)) / sd(nBGOverlaps)
  pval.perm = sum(nBGOverlaps > nLOverlaps)/length(nBGOverlaps)
  
  d <- data.frame(
    # nCADSNPs = length(o),
    Z = z_score,
    pval.z = pnorm(-abs(z_score)), ## one-tailed test ,
    signed.log10p = -log10( pnorm(-abs(z_score))) * sign(z_score), 
    pval.perm = pval.perm,
    signed.log10pperm = -log10(pval.perm) * sign(pval.perm),
    Obs.overlaps = nLOverlaps,
    OR = OR
  )
  
  if(!return_bg){return(d)}
  else{
    return(list(d,nBGOverlaps))
  }
}

plot_enrichment <- function(ztest.list, ## list object returned by Peak Z test
                            perm.p = F, ## whether to use permutation-based p values
                            binwidth = 4){
  require(ggplot2)
  if(!is.list(ztest.list)) {
    stop("Input must be a list!")
  }
  stopifnot(length(ztest.list) >1)
  if(! is.ggplot(ztest.list[[2]]) ){
    stop("Make sure a ggplot is returned by setting return_bg to TRUE.")
  }
  
  if(!perm.p){
    pval <- round(ztest.list[[1]]$pval.z,4)
  }else{
    pval <- round(ztest.list[[1]]$pval.perm,4)
  }
  if(pval == 0){
    pval.label = "p-value < 0.0001"
  }else{
    pval.label = paste("p-value: ", pval, sep = "")
  }
  nBGOverlaps <- ztest.list[[2]]
  p <- ggplot(data.frame(nOverlaps=nBGOverlaps), aes(x=nOverlaps)) +
    geom_histogram(aes(y=..density..), binwidth = binwidth) +
    geom_density()+
    geom_vline(xintercept = ztest.list[[1]]$Obs.overlaps, color = "red", linetype="dashed", size=2)+
    ggtitle(label = pval.label)+
    theme_bw(base_size = 15,base_rect_size = 2, base_line_size = 0.2) + 
    theme(plot.title = element_text(color = "red", size = 12, face = "bold", hjust = 1))
  return(p)
}
