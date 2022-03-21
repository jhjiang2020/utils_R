calc_knn_cos_sim <- function(exp, pos, k, n_null = 100, n_PCs = 20) {
  
  # calculate binary weight matrix of k-nearest neighbors
  dist_sq <- outer(rowSums(pos^2), rowSums(pos^2), "+") - 2 * tcrossprod(pos, pos)
  dist_sq[dist_sq < 0] <- 0 # numeric stability
  
  weights <- apply(dist_sq, 1, function(x) {
    thresholds <- unique(sort(x, decreasing = F))[k]
    (x <= thresholds) * 1.0 %>% return()
  }) %>% t()
  
  # calculate cosine similarity use the first n_PCs pcs
  pca <- prcomp(exp %>% t())$x[, 1:n_PCs]
  cos_sim <- tcrossprod(pca) / (sqrt(tcrossprod(rowSums(pca^2))))
  
  # extract cosine similarity between neighbors
  cos_weights <- weights * cos_sim
  cos_weights_vec <- cos_weights[weights != 0]
  # print(sum(weights != 0))
  
  # calculate background cosine similarity level
  null_weights_vec <- c()
  
  for (i in 1:n_null) {
    new_row <- sample(1:ncol(weights))
    new_col <- sample(1:ncol(weights))
    new_weights <- weights[new_row, new_col]
    
    null_weights <- new_weights * cos_sim
    null_weights_vec <- c(null_weights_vec, null_weights[new_weights != 0])
  }
  # return quantiles
  q <- seq(0, 1, 0.01)
  data.frame(
    quantile = q,
    obs = quantile(cos_weights_vec, q),
    null = quantile(null_weights_vec, q)
  ) %>%
    gather("from", "cum", -c(quantile)) %>%
    return()
}

# calculate pairwise cosine similarity per sample
calc_knn_cos_sim_libd <- function(spe, sample_id, k, n_null = 100) {
  sample_idx <- colData(spe)$sample_id == sample_id
  
  # extract spatial coordinates
  pos <- spatialData(spe)[sample_idx, c("array_row", "array_col")] %>% as.matrix()
  
  # calculate binary weight matrix of k-nearest neighbors
  dist_sq <- outer(rowSums(pos^2), rowSums(pos^2), "+") - 2 * tcrossprod(pos, pos)
  weights <- apply(dist_sq, 1, function(x) {
    thresholds <- unique(sort(x, decreasing = F))[k]
    (x <= thresholds) * 1.0 %>% return()
  }) %>% t()
  
  # calculate cosine similarity use the first 20 pcs
  pca <- reducedDim(spe)[sample_idx, 1:20]
  cos_sim <- tcrossprod(pca) / (sqrt(tcrossprod(rowSums(pca^2))))
  
  # extract cosine similarity between neighbors
  cos_weights <- weights * cos_sim
  cos_weights_vec <- cos_weights[weights != 0]
  # print(sum(weights!=0))
  
  # calculate background cosine similarity level
  null_weights_vec <- c()
  
  for (i in 1:n_null) {
    new_row <- sample(1:ncol(weights))
    new_col <- sample(1:ncol(weights))
    new_weights <- weights[new_row, new_col]
    
    null_weights <- new_weights * cos_sim
    null_weights_vec <- c(null_weights_vec, null_weights[null_weights != 0])
  }
  
  # return quantiles
  q <- seq(0, 1, 0.01)
  data.frame(
    sample = sample_id,
    quantile = q,
    obs = quantile(cos_weights_vec, q),
    null = quantile(null_weights_vec, q)
  ) %>%
    gather("from", "cum", -c(sample, quantile)) %>%
    return()
}

calc_cos_sim <- function(exp, pos, max_k = 50, topk = TRUE, n_PCs = 20) {
  
  # calculate threshold of k-nearest neighbors
  dist_sq <- outer(rowSums(pos^2), rowSums(pos^2), "+") - 2 * tcrossprod(pos, pos)
  dist_sq[dist_sq < 0] <- 0 # numeric stability
  thresholds <- apply(dist_sq, 1, function(x) {
    unique(sort(x, decreasing = F))
  })
  
  # calculate cosine similarity use the first n_PCs pcs
  pca <- prcomp(exp %>% t())$x[, 1:n_PCs]
  cos_sim <- tcrossprod(pca) / (sqrt(tcrossprod(rowSums(pca^2))))
  
  lapply(1:max_k, function(k) {
    # calculate cos_sim distribution for each k
    lapply(1:nrow(pos), function(i) {
      index <- dist_sq[i, ] <= thresholds[[i]][k]
      if ((!topk) & k > 1) {
        index <- index & (dist_sq[i, ] > thresholds[[i]][k - 1])
      }
      
      cos_sim[i, which(index)] %>% return()
    }) %>%
      do.call(what = "c") %>%
      quantile(c(0.1, 0.25, 0.5, 0.75, 0.9)) %>%
      t() %>%
      data.frame() %>%
      `names<-`(paste0("Q", c(10, 25, 50, 75, 90))) %>%
      mutate(k = k) %>%
      return()
  }) %>%
    do.call(what = "rbind") %>%
    return()
}


calc_cos_sim_libd <- function(spe, max_k = 50, topk = TRUE) {
  samples <- unique(colData(spe)$sample_id)
  
  list_sim_k <- lapply(1:max_k, function(x) NULL) %>%
    `names<-`(paste0("K_", 1:max_k))
  
  for (id in samples) {
    sample_idx <- colData(spe)$sample_id == id
    # extract spatial coordinates
    pos <- spatialData(spe)[sample_idx, c("array_row", "array_col")] %>% as.matrix()
    
    # calculate threshold of k-nearest neighbors
    dist_sq <- outer(rowSums(pos^2), rowSums(pos^2), "+") - 2 * tcrossprod(pos, pos)
    dist_sq[dist_sq < 0] <- 0 # numeric stability
    thresholds <- apply(dist_sq, 1, function(x) {
      unique(sort(x, decreasing = F))
    })
    
    # calculate cosine similarity use the first 20 pcs
    pca <- reducedDim(spe)[sample_idx, 1:20]
    cos_sim <- tcrossprod(pca) / (sqrt(tcrossprod(rowSums(pca^2))))
    
    list_sim_k_s <- lapply(1:max_k, function(k) {
      # calculate cos_sim distribution for each k
      lapply(1:nrow(pos), function(i) {
        index <- dist_sq[i, ] <= thresholds[[i]][k]
        if ((!topk) & k > 1) {
          index <- index & (dist_sq[i, ] > thresholds[[i]][k - 1])
        }
        
        cos_sim[i, which(index)] %>% return()
      }) %>%
        do.call(what = "c") %>%
        return()
    }) %>% `names<-`(paste0("K_", 1:max_k))
    
    list_sim_k <- Map(c, list_sim_k, list_sim_k_s)
  }
  
  lapply(1:max_k, function(k) {
    list_sim_k[[k]] %>%
      quantile(c(0.1, 0.25, 0.5, 0.75, 0.9)) %>%
      t() %>%
      data.frame() %>%
      `names<-`(paste0("Q", c(10, 25, 50, 75, 90))) %>%
      mutate(k = k) %>%
      return()
  }) %>%
    do.call(what = "rbind") %>%
    return()
}

estimate_decay_rate <- function(df, x, y, ignore_self = TRUE, return_model = FALSE){
  df$x = df[[x]]
  df$y = df[[y]]
  
  if (ignore_self){
    df = df %>% filter(x != 1) %>% mutate(x = x-1)
  }
  
  # estimate initial position
  y0 <- min(df$y) - 0.01
  model0 <- lm(log(y - y0) ~ x, data=df)
  alpha0 <- exp(coef(model0)[1])
  beta0 <- coef(model0)[2]
  start <- list(alpha = alpha0, beta = beta0, theta = y0)
  
  # estimate decay rate
  nlc <- nls.control(maxiter = 1000)
  model <- nls(y ~ alpha * exp(beta * x) + theta ,
               data = df, start = start,
               control = nlc)
  if (return_model){
    return(model)
  } else{
    return(coef(model)['beta'])
  }
}

plot_k_sim <- function(df_k, title = "", k = 6) {
  df_k %>%
    mutate(from = recode(from,
                         obs = "Observed",
                         null = "Shuffled"
    )) %>%
    ggplot(aes(y = cum, x = quantile, group = from, color = from)) +
    # facet_wrap(~sample) +
    geom_line() +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray") +
    labs(
      x = "Quantile",
      y = sprintf("Cumulative density of pairwise cosine similarity\nbetween neighboring spots (k=%s)", k),
      color = "",
      title = title
    ) +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    ggsci::scale_color_npg() +
    theme_classic()
}

plot_all_sim <- function(df_all, title = "", plot_type = "shadow") {
  plot_type <- match.arg(plot_type, c("shadow", "boxplot", "lineplot"))
  if (plot_type == "shadow") {
    g <- df_all %>% ggplot(aes(x = k)) +
      geom_ribbon(aes(ymin = Q10, ymax = Q90), alpha = 0.2) +
      geom_ribbon(aes(ymin = Q25, ymax = Q75), alpha = 0.3) +
      geom_line(aes(x = k, y = Q50), color = "dark red")
  } else if (plot_type == "boxplot") {
    g <- df_all %>% ggplot(aes(x = k)) +
      geom_segment(aes(x = k, xend = k, y = Q10, yend = Q90)) +
      geom_rect(
        aes(xmin = k - 0.4, xmax = k + 0.4, ymin = Q25, ymax = Q75),
        fill = "white", color = "black"
      ) +
      geom_segment(aes(x = k - 0.4, xend = k + 0.4, y = Q50, yend = Q50), color = "dark red") +
      geom_segment(aes(x = k - 0.4, xend = k + 0.4, y = Q10, yend = Q10)) +
      geom_segment(aes(x = k - 0.4, xend = k + 0.4, y = Q90, yend = Q90))
  } else {
    df_all_2 <- df_all %>% gather(key = "quantile", value = "sim", -k)
    
    g <- df_all_2 %>%
      ggplot(aes(x = k, group = quantile)) +
      geom_line(aes(y = sim, color = quantile)) +
      geom_text(
        data = df_all_glioma2 %>% filter(k == last(k)),
        aes(
          label = quantile, x = k + 0.5, y = sim + 0.05,
          color = quantile
        )
      )
  }
  
  g +
    labs(
      x = "Neighborhood degree", y = "Pairwise cosine similarity",
      title = title
    ) +
    theme_classic() +
    theme(legend.position = "none") %>% return()
}