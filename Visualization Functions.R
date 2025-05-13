### ==================================================== ###
### plotSCCpanel — 3-panel SCC brain visualization       ###
### ==================================================== ###

plotSCCpanel <- function(scc, title = NULL, zlim = "auto", palette = "nih", showAxis = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Package 'patchwork' is required.")
  if (!requireNamespace("fields", quietly = TRUE)) stop("Package 'fields' is required.")
  if (!requireNamespace("viridis", quietly = TRUE)) stop("Package 'viridis' is required.")
  
  # 1. Data prep
  Z <- scc$Z[scc$ind.inside.cover, , drop = FALSE]
  lower <- scc$scc[, 1, 1]
  upper <- scc$scc[, 1, 2]
  mean  <- scc$scc[, 1, 3]
  if (length(lower) != nrow(Z)) stop("Mismatch between SCC values and coordinates.")
  if (is.null(title)) title <- "SCC FOR GROUP 1"
  
  # 2. Axes
  xDim <- max(Z[, 1], na.rm = TRUE)
  yDim <- max(Z[, 2], na.rm = TRUE)
  xlab <- paste0("Horizontal (0-", xDim, ")")
  ylab <- paste0("Longitudinal (0-", yDim, ")")
  
  # 3. z-limits and color
  all_values <- c(lower, mean, upper)
  if (identical(zlim, "auto")) zlim <- range(all_values, na.rm = TRUE)
  palette_colors <- switch(palette,
                           "nih"     = fields::tim.colors(64),
                           "viridis" = viridis::viridis(64),
                           "gray"    = gray.colors(64),
                           stop("Invalid palette.")
  )
  fill_scale <- ggplot2::scale_fill_gradientn(
    colors = palette_colors, limits = zlim, oob = scales::squish, name = "value"
  )
  
  # 4. Plot builder
  build_tile <- function(vals, label, title_x = NULL, title_y = NULL, show_legend = FALSE) {
    df <- data.frame(x = Z[, 1], y = Z[, 2], value = vals)
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::geom_tile(show.legend = show_legend) +
      fill_scale +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = label, x = title_x, y = title_y) +
      ggplot2::theme_minimal(base_family = "serif") +
      ggplot2::theme(
        axis.title.x = if (!is.null(title_x)) ggplot2::element_text(size = 12, margin = ggplot2::margin(t = 10)) else ggplot2::element_blank(),
        axis.title.y = if (!is.null(title_y)) ggplot2::element_text(size = 12, margin = ggplot2::margin(r = 10)) else ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(size = 9),
        axis.text.y = ggplot2::element_text(size = 9),
        axis.ticks = ggplot2::element_line(size = 0.2),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 12, margin = ggplot2::margin(b = 5))
      )
  }
  
  # 5. Assemble layout
  p1 <- build_tile(lower, "Lower Band",           title_y = ylab, show_legend = FALSE)
  p2 <- build_tile(mean,  "Mean Estimate", title_x = xlab, show_legend = TRUE)
  p3 <- build_tile(upper, "Upper Band",                         show_legend = FALSE)
  
  patchwork::wrap_plots(p1, p2, p3, nrow = 1, guides = "collect") +
    patchwork::plot_annotation(
      title = toupper(title),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(family = "serif", face = "bold", size = 16, hjust = 0.5, margin = ggplot2::margin(b = 12))
      )
    )
}

### ==================================================== ###
### plotSCCcomparisonPanel — group vs group SCC panel    ###
### ==================================================== ###

plotSCCcomparisonPanel <- function(scc,
                                   title = NULL,
                                   zlim = "auto",
                                   palette = "nih",
                                   label1 = "Group 1 Mean",
                                   label2 = "Group 2 Mean",
                                   label3 = "SCC Overlay",
                                   overlay = c("positive", "negative", "both", "none"),
                                   useRawMeans = FALSE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required.")
  if (!requireNamespace("fields", quietly = TRUE)) stop("fields is required.")
  if (!requireNamespace("viridis", quietly = TRUE)) stop("viridis is required.")
  if (!requireNamespace("scales", quietly = TRUE)) stop("scales is required.")
  if (!requireNamespace("neuroSCC", quietly = TRUE)) stop("neuroSCC is required.")
  
  overlay <- match.arg(overlay)
  
  # 1. Extract coordinates and means
  if (!is.null(scc$Yhat) && !useRawMeans) {
    Z <- scc$Z[scc$ind.inside.cover, , drop = FALSE]
    mean1_vals <- scc$Yhat[1, ]
    mean2_vals <- scc$Yhat[2, ]
  } else {
    Z <- scc$Z.band
    mean1_vals <- colMeans(scc$Ya)
    mean2_vals <- colMeans(scc$Yb)
  }
  
  df1 <- data.frame(x = Z[, 1], y = Z[, 2], value = mean1_vals)
  df2 <- data.frame(x = Z[, 1], y = Z[, 2], value = mean2_vals)
  df3 <- df2  # grayscale PET background
  
  # 2. Overlay points
  overlayPts <- neuroSCC::getPoints(scc)
  posPts <- overlayPts$positivePoints
  negPts <- overlayPts$negativePoints
  
  # 3. Axis labels and title
  xDim <- max(Z[, 1], na.rm = TRUE)
  yDim <- max(Z[, 2], na.rm = TRUE)
  xlab <- paste0("Horizontal (0-", xDim, ")")
  ylab <- paste0("Longitudinal (0-", yDim, ")")
  if (is.null(title)) title <- "SCC COMPARISON PANEL"
  
  # 4. Color palette and limits
  palette_colors <- switch(palette,
                           "nih"     = fields::tim.colors(64),
                           "viridis" = viridis::viridis(64),
                           "gray"    = gray.colors(64),
                           stop("Invalid palette.")
  )
  all_values <- c(df1$value, df2$value)
  if (identical(zlim, "auto")) zlim <- range(all_values, na.rm = TRUE)
  
  fill_scale <- ggplot2::scale_fill_gradientn(
    colors = palette_colors,
    limits = zlim,
    oob = scales::squish,
    name = "value"
  )
  
  # 5. Internal tile builder
  build_tile <- function(df, label, title_x = NULL, title_y = NULL, show_legend = FALSE, scale = fill_scale) {
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::geom_tile(show.legend = show_legend) +
      scale +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = label, x = title_x, y = title_y) +
      ggplot2::theme_minimal(base_family = "serif") +
      ggplot2::theme(
        axis.title.x = if (!is.null(title_x)) ggplot2::element_text(size = 12, margin = ggplot2::margin(t = 10)) else ggplot2::element_blank(),
        axis.title.y = if (!is.null(title_y)) ggplot2::element_text(size = 12, margin = ggplot2::margin(r = 10)) else ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_text(size = 9),
        axis.text.y  = ggplot2::element_text(size = 9),
        axis.ticks   = ggplot2::element_line(size = 0.2),
        plot.title   = ggplot2::element_text(hjust = 0.5, size = 12, margin = ggplot2::margin(b = 5))
      )
  }
  
  # 6. SCC overlay over grayscale PET
  p3 <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = df3, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::scale_fill_gradient(low = "black", high = "white", name = NULL) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = label3, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 9),
      axis.text.y = ggplot2::element_text(size = 9),
      axis.ticks = ggplot2::element_line(size = 0.2),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, margin = ggplot2::margin(b = 5)),
      legend.position = "none"
    )
  
  if (overlay %in% c("positive", "both") && nrow(posPts) > 0) {
    p3 <- p3 + ggplot2::geom_point(data = posPts, ggplot2::aes(x = x, y = y), inherit.aes = FALSE,
                                   color = "blue", size = 1.6, shape = 15)
  }
  if (overlay %in% c("negative", "both") && nrow(negPts) > 0) {
    p3 <- p3 + ggplot2::geom_point(data = negPts, ggplot2::aes(x = x, y = y), inherit.aes = FALSE,
                                   color = "red", size = 1.6, shape = 17)
  }
  
  # 7. Assemble layout
  p1 <- build_tile(df1, label = label1, title_y = ylab, show_legend = FALSE)
  p2 <- build_tile(df2, label = label2, title_x = xlab, show_legend = TRUE)
  
  patchwork::wrap_plots(p1, p2, p3, nrow = 1, guides = "collect") +
    patchwork::plot_annotation(
      title = toupper(title),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(family = "serif", face = "bold", size = 16, hjust = 0.5, margin = ggplot2::margin(b = 12))
      )
    )
}


### ==================================================== ###
### plotValidationPanel — ROI vs SCC vs SPM overlay plot  ###
### ==================================================== ###

plotValidationPanel <- function(template,
                                backgroundMatrix,
                                roiPoints,
                                sccPoints,
                                spmPoints,
                                title = "Performance Validation Panel",
                                label1 = "Ground Truth (ROI)",
                                label2 = "SCC Detected",
                                label3 = "SPM Detected",
                                overlay = c("positive", "negative", "both", "none")) {
  
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  stopifnot(requireNamespace("patchwork", quietly = TRUE))
  stopifnot(requireNamespace("neuroSCC", quietly = TRUE))
  
  overlay <- match.arg(overlay)
  
  dims <- neuroSCC::getDimensions(template)
  Z <- as.matrix(expand.grid(y = 1:dims$yDim, x = 1:dims$xDim)[, c(2, 1)])
  if (ncol(backgroundMatrix) != nrow(Z)) {
    stop("backgroundMatrix must have columns equal to number of voxels in slice.")
  }
  
  dfBase <- data.frame(x = Z[, 1], y = Z[, 2], value = colMeans(backgroundMatrix))
  
  xlab <- paste0("Horizontal (0-", dims$xDim, ")")
  ylab <- paste0("Longitudinal (0-", dims$yDim, ")")
  
  fill_scale <- ggplot2::scale_fill_gradient(low = "black", high = "white", name = NULL)
  
  base_theme <- ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background  = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid       = ggplot2::element_blank(),
      axis.text        = ggplot2::element_text(size = 9),
      axis.ticks       = ggplot2::element_line(size = 0.2),
      legend.position  = "none"
    )
  
  # --- Panel 1: ROI Ground Truth ---
  p1 <- ggplot2::ggplot(dfBase, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_tile() +
    fill_scale +
    ggplot2::geom_point(data = roiPoints, ggplot2::aes(x = x, y = y),
                        inherit.aes = FALSE, color = "green", size = 1.2, shape = 15) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = label1, x = NULL, y = ylab) +
    base_theme +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
      plot.title.position = "panel",
      axis.title.y = ggplot2::element_text(size = 12, margin = ggplot2::margin(r = 10))
    )
  
  # --- Panel 2: SCC Detected ---
  p2 <- ggplot2::ggplot(dfBase, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_tile() +
    fill_scale +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = label2, x = xlab, y = NULL) +
    base_theme +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
      plot.title.position = "panel",
      axis.title.x = ggplot2::element_text(size = 12, margin = ggplot2::margin(t = 10)),
      axis.title.y = ggplot2::element_blank()
    )
  
  if (overlay %in% c("positive", "both") && nrow(sccPoints$positivePoints) > 0) {
    p2 <- p2 + ggplot2::geom_point(data = sccPoints$positivePoints, ggplot2::aes(x = x, y = y),
                                   inherit.aes = FALSE, color = "blue", size = 1.2, shape = 15)
  }
  if (overlay %in% c("negative", "both") && nrow(sccPoints$negativePoints) > 0) {
    p2 <- p2 + ggplot2::geom_point(data = sccPoints$negativePoints, ggplot2::aes(x = x, y = y),
                                   inherit.aes = FALSE, color = "red", size = 1.2, shape = 17)
  }
  
  # --- Panel 3: SPM Detected ---
  p3 <- ggplot2::ggplot(dfBase, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_tile() +
    fill_scale +
    ggplot2::geom_point(data = spmPoints, ggplot2::aes(x = x, y = y),
                        inherit.aes = FALSE, color = "blue", size = 1.2, shape = 15) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = label3, x = NULL, y = NULL) +
    base_theme +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
      plot.title.position = "panel",
      axis.title.y = ggplot2::element_blank()
    )
  
  # --- Assemble all three panels ---
  patchwork::wrap_plots(p1, p2, p3, nrow = 1) +
    patchwork::plot_annotation(
      title = toupper(title),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(family = "serif", face = "bold", size = 16,
                                           hjust = 0.5, margin = ggplot2::margin(b = 12))
      )
    )
}
