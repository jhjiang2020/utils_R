require(ggplot2)
require(Seurat)
vega_10_palette <- c("#1F78B4","#F57E20","#2AA148","#D72C27","#9268AC","#8B554D", "#D979B1", "#7E8180", "#BCBC32","#1DBFCF")

vega_20_palette <- c("#1F78B4", "#AEC6E7", "#F57E20","#FCBA78", "#2AA148","#9FD18A", "#D72C27","#F69797", "#9268AC","#C5B0D5" ,"#8B554D","#C49C94" , "#D979B1","#F5B5D1" , "#7E8180","#C8C7C6" , "#BCBC32","#DBDB8C" ,"#1DBFCF", "#A0DBE5")

FeaturePlot_autoscale <- function(gene, object, reduction){
  max <- GetAssayData(object, assay="RNA", slot="data")[gene,] %>% max
  p <- FeaturePlot(object, reduction = reduction, features = gene, slot = "data") +
    scale_color_gradientn(limits=c(0,max),colors=c("grey95", BuenColors::jdb_palette("solar_rojos", type="continuous")), labels = c("Min", "Max"), breaks=c(0,max),
                          name = "Expression level",
                          guide = guide_colorbar(frame.linewidth = 0.5, frame.colour = "black"))
  return(p)
}