## modified from the source code of UpSetR package v1.4.0 to enable customized size bar to the left of the matrix
## environment(upset) <- asNamespace('UpSetR')
## assignInNamespace(x="upset", value=upset, ns="UpSetR")
upset <- function(data, nsets = 5, set_size_data_manual=NULL, nintersects = 40, sets = NULL, keep.order = F, set.metadata = NULL, intersections = NULL,
                  matrix.color = "gray23", main.bar.color = "gray23", mainbar.y.label = "Intersection Size", mainbar.y.max = NULL,
                  sets.bar.color = "gray23", plot.title = NA, sets.x.label = "Set Size", point.size = 2.2, line.size = 0.7,
                  mb.ratio = c(0.70,0.30), expression = NULL, att.pos = NULL, att.color = main.bar.color, order.by = c("freq", "degree"),
                  decreasing = c(T, F), show.numbers = "yes", number.angles = 0, number.colors=NULL, group.by = "degree",cutoff = NULL,
                  queries = NULL, query.legend = "none", shade.color = "gray88", shade.alpha = 0.25, matrix.dot.alpha =0.5,
                  empty.intersections = NULL, color.pal = 1, boxplot.summary = NULL, attribute.plots = NULL, scale.intersections = "identity",
                  scale.sets = "identity", text.scale = 1, set_size.angles = 0 , set_size.show = FALSE, set_size.numbers_size = NULL, set_size.scale_max = NULL){
  
  startend <-FindStartEnd(data)
  first.col <- startend[1]
  last.col <- startend[2]
  
  if(color.pal == 1){
    palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2",
                 "#7F7F7F", "#BCBD22", "#17BECF")
  }
  else{
    palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
                 "#CC79A7")
  }
  
  if(is.null(intersections) == F){
    Set_names <- unique((unlist(intersections)))
    Sets_to_remove <- Remove(data, first.col, last.col, Set_names)
    New_data <- Wanted(data, Sets_to_remove)
    Num_of_set <- Number_of_sets(Set_names)
    if(keep.order == F){
      Set_names <- order_sets(New_data, Set_names)
    }
    All_Freqs <- specific_intersections(data, first.col, last.col, intersections, order.by, group.by, decreasing,
                                        cutoff, main.bar.color, Set_names)
  }
  else if(is.null(intersections) == T){
    Set_names <- sets
    if(is.null(Set_names) == T || length(Set_names) == 0 ){
      Set_names <- FindMostFreq(data, first.col, last.col, nsets)
    }
    Sets_to_remove <- Remove(data, first.col, last.col, Set_names)
    New_data <- Wanted(data, Sets_to_remove)
    Num_of_set <- Number_of_sets(Set_names)
    if(keep.order == F){
      Set_names <- order_sets(New_data, Set_names)
    }
    All_Freqs <- Counter(New_data, Num_of_set, first.col, Set_names, nintersects, main.bar.color,
                         order.by, group.by, cutoff, empty.intersections, decreasing)
  }
  Matrix_setup <- Create_matrix(All_Freqs)
  labels <- Make_labels(Matrix_setup)
  #Chose NA to represent NULL case as result of NA being inserted when at least one contained both x and y
  #i.e. if one custom plot had both x and y, and others had only x, the y's for the other plots were NA
  #if I decided to make the NULL case (all x and no y, or vice versa), there would have been alot more if/else statements
  #NA can be indexed so that we still get the non NA y aesthetics on correct plot. NULL cant be indexed.
  att.x <- c(); att.y <- c();
  if(is.null(attribute.plots) == F){
    for(i in seq_along(attribute.plots$plots)){
      if(length(attribute.plots$plots[[i]]$x) != 0){
        att.x[i] <- attribute.plots$plots[[i]]$x
      }
      else if(length(attribute.plots$plots[[i]]$x) == 0){
        att.x[i] <- NA
      }
      if(length(attribute.plots$plots[[i]]$y) != 0){
        att.y[i] <- attribute.plots$plots[[i]]$y
      }
      else if(length(attribute.plots$plots[[i]]$y) == 0){
        att.y[i] <- NA
      }
    }
  }
  
  BoxPlots <- NULL
  if(is.null(boxplot.summary) == F){
    BoxData <- IntersectionBoxPlot(All_Freqs, New_data, first.col, Set_names)
    BoxPlots <- list()
    for(i in seq_along(boxplot.summary)){
      BoxPlots[[i]] <- BoxPlotsPlot(BoxData, boxplot.summary[i], att.color)
    }
  }
  
  customAttDat <- NULL
  customQBar <- NULL
  Intersection <- NULL
  Element <- NULL
  legend <- NULL
  EBar_data <- NULL
  if(is.null(queries) == F){
    custom.queries <- SeperateQueries(queries, 2, palette)
    customDat <- customQueries(New_data, custom.queries, Set_names)
    legend <- GuideGenerator(queries, palette)
    legend <- Make_legend(legend)
    if(is.null(att.x) == F && is.null(customDat) == F){
      customAttDat <- CustomAttData(customDat, Set_names)
    }
    customQBar <- customQueriesBar(customDat, Set_names, All_Freqs, custom.queries)
  }
  if(is.null(queries) == F){
    Intersection <- SeperateQueries(queries, 1, palette)
    Matrix_col <- intersects(QuerieInterData, Intersection, New_data, first.col, Num_of_set,
                             All_Freqs, expression, Set_names, palette)
    Element <- SeperateQueries(queries, 1, palette)
    EBar_data <-ElemBarDat(Element, New_data, first.col, expression, Set_names,palette, All_Freqs)
  }
  else{
    Matrix_col <- NULL
  }
  
  Matrix_layout <- Create_layout(Matrix_setup, matrix.color, Matrix_col, matrix.dot.alpha)
  if(is.null(set_size_data_manual)){
    Set_sizes <- FindSetFreqs(New_data, first.col, Num_of_set, Set_names, keep.order)
  }else{
    Set_sizes <- set_size_data_manual
  }
  Bar_Q <- NULL
  if(is.null(queries) == F){
    Bar_Q <- intersects(QuerieInterBar, Intersection, New_data, first.col, Num_of_set, All_Freqs, expression, Set_names, palette)
  }
  QInter_att_data <- NULL
  QElem_att_data <- NULL
  if((is.null(queries) == F) & (is.null(att.x) == F)){
    QInter_att_data <- intersects(QuerieInterAtt, Intersection, New_data, first.col, Num_of_set, att.x, att.y,
                                  expression, Set_names, palette)
    QElem_att_data <- elements(QuerieElemAtt, Element, New_data, first.col, expression, Set_names, att.x, att.y,
                               palette)
  }
  AllQueryData <- combineQueriesData(QInter_att_data, QElem_att_data, customAttDat, att.x, att.y)
  
  ShadingData <- NULL
  
  if(is.null(set.metadata) == F){
    ShadingData <- get_shade_groups(set.metadata, Set_names, Matrix_layout, shade.alpha)
    output <- Make_set_metadata_plot(set.metadata, Set_names)
    set.metadata.plots <- output[[1]]
    set.metadata <- output[[2]]
    
    if(is.null(ShadingData) == FALSE){
      shade.alpha <- unique(ShadingData$alpha)
    }
  } else {
    set.metadata.plots <- NULL
  }
  if(is.null(ShadingData) == TRUE){
    ShadingData <- MakeShading(Matrix_layout, shade.color)
  }
  Main_bar <- suppressMessages(Make_main_bar(All_Freqs, Bar_Q, show.numbers, mb.ratio, customQBar, number.angles, number.colors, EBar_data, mainbar.y.label,
                                             mainbar.y.max, scale.intersections, text.scale, attribute.plots, plot.title))
  Matrix <- Make_matrix_plot(Matrix_layout, Set_sizes, All_Freqs, point.size, line.size,
                             text.scale, labels, ShadingData, shade.alpha)
  Sizes <- Make_size_plot(Set_sizes, sets.bar.color, mb.ratio, sets.x.label, scale.sets, text.scale, set_size.angles,set_size.show,
                          set_size.scale_max, set_size.numbers_size)
  
  # Make_base_plot(Main_bar, Matrix, Sizes, labels, mb.ratio, att.x, att.y, New_data,
  #                expression, att.pos, first.col, att.color, AllQueryData, attribute.plots,
  #                legend, query.legend, BoxPlots, Set_names, set.metadata, set.metadata.plots)
  
  structure(class = "upset",
            .Data=list(
              Main_bar = Main_bar,
              Matrix = Matrix,
              Sizes = Sizes,
              labels = labels,
              mb.ratio = mb.ratio,
              att.x = att.x,
              att.y = att.y,
              New_data = New_data,
              expression = expression,
              att.pos = att.pos,
              first.col = first.col,
              att.color = att.color,
              AllQueryData = AllQueryData,
              attribute.plots = attribute.plots,
              legend = legend,
              query.legend = query.legend,
              BoxPlots = BoxPlots,
              Set_names = Set_names,
              set.metadata = set.metadata,
              set.metadata.plots = set.metadata.plots)
  )
}

