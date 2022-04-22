require(patchwork)
require(GenomicRanges)
require(ggplot2)
require(ggrepel)

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
  SNP.intersect <- subsetByOverlaps(x = SNP, ranges = region)
  SNP.df <- as.data.frame(x = SNP.intersect)
  SNP.df$rsid <- label
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
        mapping = aes(x = (start+end)/2, y = -0.5, label = rsid),
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

