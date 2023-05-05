require(patchwork)
require(GenomicRanges)
require(ggplot2)
require(ggrepel)
library(fastmatch)

SNPPlot <- function(
  region,
  SNP,
  label,
  width = NULL,
  color = "dimgrey"
) {
  
  stopifnot(length(label) == length(SNP))
  
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  
  width <- round(as.numeric(width %||% width(region)/300))
  
  ## convert SNP from a character vector to GRanges vector
  if (!inherits(x = SNP, what = "GRanges")) {
    SNP <- StringToGRanges(regions = SNP) + width
  }
  
  
  # subset to covered range
  SNP$rsid <- label
  SNP.df <- as.data.frame(subsetByOverlaps(x = SNP, ranges = region))
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  
  if (nrow(x = SNP.df) > 0) {
    
    SNP.df$start[SNP.df$start < start.pos] <- start.pos
    SNP.df$end[SNP.df$end > end.pos] <- end.pos
    SNP.plot <- ggplot(
      data = SNP.df,
      #aes_string(color = SetIfNull(x = group.by, y = "color"))
    ) +
      geom_segment(aes(x = (start+end)/2, y = -0.3, xend = (start+end)/2, yend = 0.3, color=color),
                   size = 1,
                   data = SNP.df) +
      geom_text_repel(
        data = SNP.df,
        mapping = aes(x = (start+end)/2, y = 0, label = rsid),
        min.segment.length = 0,
        position = position_dodge( width = 7),
        size = 3
      )+
      #scale_y_discrete()
      scale_y_continuous(limits = c(-1,1))
  } else {
    # no SNPs present in region, make empty panel
    SNP.plot <- ggplot(data = SNP.df)
  }
  SNP.plot <- SNP.plot + theme_classic() +
    ylab(label = "SNPs") +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(c(start.pos, end.pos)) +
    NoLegend()

  return(SNP.plot)
}

CombineTracks.2 <- function(
  plotlist,
  expression.plot = NULL,
  heights = NULL,
  widths = NULL
) {
  # remove any that are NULL
  nullplots <- sapply(X = plotlist, FUN = is.null)
  plotlist <- plotlist[!nullplots]
  heights <- heights[!nullplots]
  
  if (length(x = plotlist) == 1) {
    return(plotlist[[1]])
  }
  
  
  # remove x-axis from all but last plot
  for (i in 1:(length(x = plotlist) - 1)) {
    plotlist[[i]] <- plotlist[[i]] + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x.bottom = element_blank(),
      axis.ticks.x.bottom = element_blank()
    ) 
  }
  # combine plots
  if (is.null(x = heights)) {
    # set height of the last element to 8x more than other elements
    n.plots <- length(x = plotlist)
    heights <- c(rep(1, n.plots - 1), 8)
  } else {
    if (length(x = heights) != length(x = plotlist)) {
      stop("Relative height must be supplied for each plot")
    }
  }
  if (!is.null(x = expression.plot)) {
    for (i in 1:(length(x = plotlist) - 1)) {
      plotlist[[i]] <- plotlist[[i]] + theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x.bottom = element_blank(),
        axis.ticks.x.bottom = element_blank()
      ) + guide_area() + plot_layout(widths = widths)
    }

    # align expression plot with the last element in plot list
    n <- length(x = plotlist)
    plotlist[[n]] <-  plotlist[[n]] + theme_bw(base_size = 15,base_rect_size = 2, base_line_size = 0.2) + NoLegend()
    expression.plot <- expression.plot + theme_bw(base_size = 15,base_rect_size = 2, base_line_size = 0.2) + NoLegend()
    plotlist[[n]] <- (plotlist[[n]] + expression.plot)+ 
      plot_layout(widths = widths)
    
    p <- wrap_plots(plotlist, ncol = 1, heights = heights)
    

  } else {
    p <- wrap_plots(plotlist, ncol = 1, heights = heights)
  }
  return(p)
}

LookupGeneCoords <- function(object, gene, assay = NULL) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay")
  }
  annotations <- Annotation(object = object[[assay]])
  isgene <- annotations$gene_name == gene
  isgene <- !is.na(x = isgene) & isgene
  annot.sub <- annotations[isgene]
  if (length(x = annot.sub) == 0) {
    return(NULL)
  } else {
    gr <- GRanges(seqnames = as.character(x = seqnames(x = annot.sub))[[1]],
                  ranges = IRanges(start = min(start(x = annot.sub)),
                                   end = max(end(x = annot.sub))))
    return(gr)
  }
}

SetIfNull <- function(x, y) {
  if (is.null(x = x)) {
    return(y)
  } else {
    return(x)
  }
}

PeakPlot.2 <- function(
  region,
  peaks, # a GRanges object containing peaks and group.by information
  group.by = NULL,
  color = NULL
) {
  # subset to covered range
  peak.intersect <- subsetByOverlaps(x = peaks, ranges = region)
  peak.df <- as.data.frame(x = peak.intersect)
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  
  if (nrow(x = peak.df) > 0) {
    if (!is.null(x = group.by)) {
      if (!(group.by %in% colnames(x = peak.df))) {
        warning("Requested grouping variable not found")
        group.by <- NULL
      }
    }
    peak.df$start[peak.df$start < start.pos] <- start.pos
    peak.df$end[peak.df$end > end.pos] <- end.pos
    peak.plot <- ggplot(
      data = peak.df,
      aes_string(color = SetIfNull(x = group.by, y = "color"))
    ) +
      geom_segment(aes(x = start, y = 0, xend = end, yend = 0),
                   size = 2,
                   data = peak.df)
  } else {
    # no peaks present in region, make empty panel
    peak.plot <- ggplot(data = peak.df)
  }
  peak.plot <- peak.plot + theme_classic() +
    ylab(label = "Peaks") +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(c(start.pos, end.pos))
  if (!is.null(x = group.by)) {
    # remove legend, change color
    peak.plot <- peak.plot +
      scale_color_manual(values = color) +
      theme(legend.position = "none")
  }
  return(peak.plot)
}

AnnotationPlot.2 <- function(object, region, mode = "gene") {
  if(mode == "gene") {
    collapse_transcript <- TRUE
    label <- "gene_name"
  } else if (mode == "transcript") {
    collapse_transcript <- FALSE
    label <- "tx_id"
  } else {
    stop("Unknown mode requested, choose either 'gene' or 'transcript'")
  }
  
  if(!inherits(x = object, what = "Seurat")){
    annotation <- object
  }else{
    annotation <- Annotation(object = object)
  }
  if (is.null(x = annotation)) {
    return(NULL)
  }
  
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  
  # get names of genes that overlap region, then subset to include only those
  # genes. This avoids truncating the gene if it runs outside the region
  annotation.subset <- subsetByOverlaps(x = annotation, ranges = region)
  if (mode == "gene") {
    genes.keep <- unique(x = annotation.subset$gene_name)
    annotation.subset <- annotation[
      fmatch(x = annotation$gene_name, table = genes.keep, nomatch = 0L) > 0L
    ]
  } else {
    tx.keep <- unique(x = annotation.subset$tx_id)
    annotation.subset <- annotation[
      fmatch(x = annotation$tx_id, table = tx.keep, nomatch = 0L) > 0L
    ]
  }
  
  if (length(x = annotation.subset) == 0) {
    # make empty plot
    p <- ggplot(data = data.frame())
    y_limit <- c(0, 1)
  } else {
    annotation_df_list <- reformat_annotations(
      annotation = annotation.subset,
      start.pos = start.pos,
      end.pos = end.pos,
      collapse_transcript = collapse_transcript
    )
    p <- ggplot() +
      # exons
      geom_segment(
        data = annotation_df_list$exons,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$exons$dodge,
          xend = "end",
          yend = annotation_df_list$exons$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 3
      ) +
      # gene body
      geom_segment(
        data = annotation_df_list$labels,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$labels$dodge,
          xend = "end",
          yend = annotation_df_list$labels$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 1/2
      )
    if (nrow(x = annotation_df_list$plus) > 0) {
      # forward strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$plus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$plus$dodge,
          xend = "end",
          yend = annotation_df_list$plus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "last",
          type = "open",
          angle = 45,
          length = unit(x = 0.04, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/2
      )
    }
    if (nrow(x = annotation_df_list$minus) > 0) {
      # reverse strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$minus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$minus$dodge,
          xend = "end",
          yend = annotation_df_list$minus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "first",
          type = "open",
          angle = 45,
          length = unit(x = 0.04, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/2
      )
    }
    # label genes
    n_stack <- max(annotation_df_list$labels$dodge)
    annotation_df_list$labels$dodge <- annotation_df_list$labels$dodge + 0.2
    p <- p + geom_text(
      data = annotation_df_list$labels,
      mapping = aes_string(x = "position", y = "dodge", label = label),
      size = 2.5
    )
    y_limit <- c(0.9, n_stack + 0.4)
  }
  p <- p +
    theme_classic() +
    ylab("Genes") +
    xlab(label = paste0(chromosome, " position (bp)")) +
    coord_cartesian(xlim=c(start.pos, end.pos),
                    ylim=y_limit ) + 
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    scale_color_manual(values = c("darkblue", "darkgreen"))
  return(p)
}
reformat_annotations <- function(
    annotation,
    start.pos,
    end.pos,
    collapse_transcript = TRUE
) {
  total.width <- end.pos - start.pos
  tick.freq <- total.width / 50
  annotation <- annotation[annotation$type == "exon"]
  exons <- as.data.frame(x = annotation, row.names = NULL)
  if (collapse_transcript) {
    annotation <- split(
      x = annotation,
      f = annotation$gene_name
    )
  } else {
    annotation <- split(
      x = annotation,
      f = annotation$tx_id
    )
  }
  annotation <- lapply(X = annotation, FUN = as.data.frame, row.names=NULL)
  
  # add gene total start / end
  gene_bodies <- list()
  for (i in seq_along(annotation)) {
    df <- data.frame(
      seqnames = annotation[[i]]$seqnames[[1]],
      start = min(annotation[[i]]$start),
      end = max(annotation[[i]]$end),
      strand = annotation[[i]]$strand[[1]],
      tx_id = annotation[[i]]$tx_id[[1]],
      gene_name = annotation[[i]]$gene_name[[1]],
      gene_biotype = annotation[[i]]$gene_biotype[[1]],
      type = "body"
    )
    # trim any that extend beyond region
    df$start <- ifelse(
      test = df$start < start.pos,
      yes = start.pos,
      no = df$start
    )
    df$end <- ifelse(
      test = df$end > end.pos,
      yes = end.pos,
      no = df$end
    )
    breaks <- split_body(df = df, width = tick.freq)
    df <- rbind(df, breaks)
    gene_bodies[[i]] <- df
  }
  gene_bodies <- do.call(what = rbind, args = gene_bodies)
  
  # record if genes overlap
  overlap_idx <- record_overlapping(
    annotation = gene_bodies,
    min.gapwidth = 1000,
    collapse_transcript = collapse_transcript
  )
  # overlap_idx <- overlap_idx
  if (collapse_transcript) {
    gene_bodies$dodge <- overlap_idx[gene_bodies$gene_name]
    exons$dodge <- overlap_idx[exons$gene_name]
  } else {
    gene_bodies$dodge <- overlap_idx[gene_bodies$tx_id]
    exons$dodge <- overlap_idx[exons$tx_id]
  }
  
  label_df <- gene_bodies[gene_bodies$type == "body", ]
  label_df$width <- label_df$end - label_df$start
  label_df$position <- label_df$start + (label_df$width / 2)
  
  onplus <- gene_bodies[gene_bodies$strand %in% c("*", "+"), ]
  onminus <- gene_bodies[gene_bodies$strand == "-", ]
  
  return(
    list(
      "labels" = label_df,
      "exons" = exons,
      "plus" = onplus,
      "minus" = onminus
    )
  )
}

split_body <- function(df, width = 1000) {
  wd <- df$end - df$start
  nbreak <- wd / width
  if (nbreak > 1) {
    steps <- 0:(nbreak)
    starts <- (width * steps) + df$start
    starts[starts > df$end] <- NULL
  } else {
    starts <- df$end
  }
  breaks <- data.frame(
    seqnames = df$seqnames[[1]],
    start = starts,
    end = starts + 1,
    strand = df$strand[[1]],
    tx_id = df$tx_id[[1]],
    gene_name = df$gene_name[[1]],
    gene_biotype = df$gene_biotype[[1]],
    type = "arrow"
  )
  return(breaks)
}
record_overlapping <- function(
    annotation,
    min.gapwidth = 1000,
    collapse_transcript = TRUE
) {
  # convert back to granges
  annotation$strand <- "*"
  gr <- makeGRangesFromDataFrame(
    df = annotation[annotation$type == "body", ], keep.extra.columns = TRUE
  )
  # work out which ranges overlap
  collapsed <- reduce(
    x = gr, with.revmap = TRUE, min.gapwidth = min.gapwidth
  )$revmap
  idx <- seq_along(gr)
  for (i in seq_along(collapsed)) {
    mrg <- collapsed[[i]]
    for (j in seq_along(mrg)) {
      idx[[mrg[[j]]]] <- j
    }
  }
  if (collapse_transcript) {
    names(x = idx) <- gr$gene_name
  } else {
    names(x = idx) <- gr$tx_id
  }
  return(idx)
}

