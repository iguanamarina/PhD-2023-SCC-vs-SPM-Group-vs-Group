############################## ################### ############## ### 
##
## Script name: MASTERSCRIPT (for practical case)
##
## Purpose of script: Script to import NIFTI files to R, create a database, get PPT names,
## create triangulation parameters, and carry on with SCCs calculation. This script asumes
## that Matlab NIFTI pre-processing has already been carried out correctly.
##
## Date Created: 2022-01-10
## Date Updated: 2025-04-24
##
## Author: Juan A. Arias (M.Sc.)
## Email: juanantonio.arias.lopez@usc.es
## Webpage: https://juan-arias.xyz
##   
############################## ################### ############## ### 


### ==================================================== ###
###                PREAMBLE: Package Setup               ###
### ==================================================== ###

# Set project working directory
setwd("~/GitHub/PhD-2023-SCC-vs-SPM-Group-vs-Group")

# -------------------------------
# CRAN Packages
# -------------------------------
cran_packages <- c(
  # Core dependencies
  "mgcv", "gamair", "oro.nifti", "memisc",
  # neuroSCC package
  "neuroSCC",
  # Visualization
  "ggplot2", "patchwork", "fields", "viridis",
  # Tables and reporting
  "knitr", "kableExtra", "tibble", "magrittr",
  # Data wrangling
  "tidyr", "dplyr", "scales",
  # SCC-related and additional tools
  "devtools", "remotes", "readr", "imager", "itsadug", "contoureR"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}
rm(pkg, cran_packages)

# -------------------------------
# Non-CRAN Packages via drat
# -------------------------------
non_cran_packages <- c("Triangulation", "ImageSCC", "BPST")

for (pkg in non_cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(
      pkg,
      repos = c(neuroSCC = "https://iguanamarina.github.io/neuroSCC-drat",
                CRAN = "https://cloud.r-project.org"),
      type = "source"
    )
  }
  library(pkg, character.only = TRUE)
}
rm(pkg, non_cran_packages)

### ==================================================== ###
###        Testing `neuroSCC::neuroCleaner()` function   ###
### ==================================================== ###

# Load demographics table (used for optional metadata attachment)
demo <- read.csv2("PETimg_masked for practical case/Demographics.csv", stringsAsFactors = FALSE)

# Run a simple test: convert one NIfTI file into a tidy dataframe
example <- neuroSCC::neuroCleaner("PETimg_masked for practical case/003_S_1059.img", 
                                  demo = demo, demoRow = 51)

### ==================================================== ###
###    Build CN and AD Databases Based on Demographics   ###
### ==================================================== ###

# Define output paths (relative to root)
dbPathCN <- "Auxiliary Files/databaseCN.RData"
dbPathAD <- "Auxiliary Files/databaseAD.RData"

# Absolute paths for file.exists()
fullPathCN <- file.path("Auxiliary Files", "databaseCN.RData")
fullPathAD <- file.path("Auxiliary Files", "databaseAD.RData")

# Go to the folder with image files
setwd("PETimg_masked for practical case")

# List available .img files
allFiles <- list.files(pattern = "\\.img$", full.names = FALSE)

# Add .img suffix to demo entries
demo$filename <- paste0(demo$PPT, ".img")

# ---- CONTROL GROUP ----
if (!file.exists(file.path("..", fullPathCN))) {
  
  message("[INFO] Creating Control group database...")
  
  filesCN <- intersect(allFiles, demo$filename[demo$Group == "CN"])
  # demoCN <- demo[demo$filename %in% filesCN, ]
  
  databaseCN <- neuroSCC::databaseCreator(
    pattern = paste(filesCN, collapse = "|"),  # regex match
    control = TRUE,
    demo = NULL,
    quiet = FALSE
  )
  
  databaseCN$pet[databaseCN$pet <= 0] <- NaN
  
  setwd("..")
  save(databaseCN, file = fullPathCN)
  message("[INFO] Control group database saved to ", fullPathCN)
  
} else {
  message("[INFO] Control group database already exists — loading.")
  load(fullPathCN)
}

# ---- AD GROUP ----
setwd("PETimg_masked for practical case")

if (!file.exists(file.path("..", fullPathAD))) {
  
  message("[INFO] Creating AD group database...")
  
  filesAD <- intersect(allFiles, demo$filename[demo$Group == "AD"])
  demoAD <- demo[demo$filename %in% filesAD, ]
  
  databaseAD <- neuroSCC::databaseCreator(
    pattern = paste(filesAD, collapse = "|"),
    control = FALSE,
    demo = demoAD,
    quiet = FALSE
  )
  
  databaseAD$pet[databaseAD$pet <= 0] <- NaN
  
  setwd("..")
  save(databaseAD, file = fullPathAD)
  message("[INFO] AD group database saved to ", fullPathAD)
  
} else {
  message("[INFO] AD group database already exists — loading.")
  load(fullPathAD)
}

### ==================================================== ###
###       Create SCC Matrices for CN and AD Groups       ###
### ==================================================== ###

# Define Z level
paramZ <- 35

# Ensure CN and AD databases are available
if (!exists("databaseCN")) {
  message("[INFO] Loading Control group database...")
  load("Auxiliary Files/databaseCN.RData")
}

if (!exists("databaseAD")) {
  message("[INFO] Loading AD group database...")
  load("Auxiliary Files/databaseAD.RData")
}

# Define matrix file paths
matrixPathCN <- "Auxiliary Files/matrixControls.RDS"
matrixPathAD <- "Auxiliary Files/matrixAD.RDS"

# ---- CONTROL GROUP MATRIX ----
if (!file.exists(matrixPathCN)) {
  
  message("[INFO] Creating matrix for Control group (z = ", paramZ, ")...")
  
  matrixControls <- neuroSCC::matrixCreator(
    database = databaseCN,
    paramZ = paramZ,
    quiet = FALSE
  )
  
  saveRDS(matrixControls, file = matrixPathCN)
  message("[INFO] Control group matrix saved to ", matrixPathCN)
  
} else {
  message("[INFO] Control group matrix already exists — loading.")
  matrixControls <- readRDS(matrixPathCN)
}

# ---- AD GROUP MATRIX ----
if (!file.exists(matrixPathAD)) {
  
  message("[INFO] Creating matrix for AD group (z = ", paramZ, ")...")
  
  matrixAD <- neuroSCC::matrixCreator(
    database = databaseAD,
    paramZ = paramZ,
    quiet = FALSE
  )
  
  saveRDS(matrixAD, file = matrixPathAD)
  message("[INFO] AD group matrix saved to ", matrixPathAD)
  
} else {
  message("[INFO] AD group matrix already exists — loading.")
  matrixAD <- readRDS(matrixPathAD)
}


### ==================================================== ###
###     Mean Normalization for Control and AD Matrices   ###
### ==================================================== ###

message("[INFO] Applying mean normalization to Control group matrix...")
matrixControlsNormalized <- neuroSCC::meanNormalization(
  matrixData = matrixControls,
  quiet = FALSE
)

message("[INFO] Applying mean normalization to AD group matrix...")
matrixADNormalized <- neuroSCC::meanNormalization(
  matrixData = matrixAD,
  quiet = FALSE
)

### ==================================================== ###
###         Generate Triangulation Grid (n = 15)          ###
### ==================================================== ###

# Create triangulation using outer boundary only, no holes
triMesh <- Triangulation::TriMesh(contourLines[[1]], n = 15)

# Extract estimated and banding components for SCC
V.est <- as.matrix(triMesh[[1]])
Tr.est <- as.matrix(triMesh[[2]])
V.band <- V.est
Tr.band <- Tr.est

### ==================================================== ###
###      Define Parameters for SCC Estimation (Wang et al.) ###
### ==================================================== ###

d.est <- 5                           # Spline degree for mean function
d.band <- 2                          # Spline degree for SCC
r <- 1                               # Smoothing parameter
lambda <- 10^seq(-6, 3, 0.5)         # Grid of penalty parameters
alpha.grid <- c(0.10, 0.05, 0.01)    # Confidence levels

### ==================================================== ###
###     SCC Construction and Save for One-Group Models    ###
### ==================================================== ###

# Paths to precomputed result files
sccCNPath <- "z35/results/SCC_CN.RDS"
sccADPath <- "z35/results/SCC_AD.RDS"

# Create Z coordinate grid from the mask file
template <- "Auxiliary Files/new_mask.nii"
dims <- neuroSCC::getDimensions(template)
Z <- as.matrix(expand.grid(y = 1:dims$yDim, x = 1:dims$xDim)[, c(2, 1)])

# Ensure we’re in project root before saving
setwd("~/GitHub/PhD-2023-SCC-vs-SPM-Group-vs-Group")

# ---- Control Group ----
if (!file.exists(sccCNPath)) {
  message("[INFO] Running SCC estimation for Control group...")
  
  sccCN <- ImageSCC::scc.image(
    Ya = matrixControlsNormalized,
    Z = Z,
    d.est = d.est,
    d.band = d.band,
    r = r,
    V.est.a = V.est, Tr.est.a = Tr.est,
    V.band.a = V.band, Tr.band.a = Tr.band,
    penalty = TRUE,
    lambda = lambda,
    alpha.grid = alpha.grid,
    adjust.sigma = TRUE
  )
  
  saveRDS(sccCN, file = sccCNPath)
  message("[INFO] Control group SCC result saved to ", sccCNPath)
} else {
  message("[INFO] Control group SCC result already exists — skipping.")
}

# ---- AD Group ----
if (!file.exists(sccADPath)) {
  message("[INFO] Running SCC estimation for AD group...")
  
  sccAD <- ImageSCC::scc.image(
    Ya = matrixADNormalized,
    Z = Z,
    d.est = d.est,
    d.band = d.band,
    r = r,
    V.est.a = V.est, Tr.est.a = Tr.est,
    V.band.a = V.band, Tr.band.a = Tr.band,
    penalty = TRUE,
    lambda = lambda,
    alpha.grid = alpha.grid,
    adjust.sigma = TRUE
  )
  
  saveRDS(sccAD, file = sccADPath)
  message("[INFO] AD group SCC result saved to ", sccADPath)
} else {
  message("[INFO] AD group SCC result already exists — skipping.")
}

# Custom function plotSCCpanel()

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

# Plot SCC results using the custom function
plotSCCpanel(sccCN, title = "Control Subject")
plotSCCpanel(sccAD, title = "AD Subject")


### ==================================================== ###
###      SCC Comparison Between AD and CN Groups (z = 35) ###
### ==================================================== ###

# Define save path
sccComparisonPath <- "z35/results/SCC_Comparison_AD_vs_CN.RDS"

# Ensure project root working directory
setwd("~/GitHub/PhD-2023-SCC-vs-SPM-Group-vs-Group")

# Compute coordinate matrix Z
template <- "Auxiliary Files/new_mask.nii"
dims <- neuroSCC::getDimensions(template)
Z <- as.matrix(expand.grid(y = 1:dims$yDim, x = 1:dims$xDim)[, c(2, 1)])

# ---- SCC Group Comparison ----
if (!file.exists(sccComparisonPath)) {
  
  message("[INFO] Running SCC group comparison (AD vs CN)...")
  
  sccGroupComparison <- ImageSCC::scc.image(
    Ya = matrixADNormalized,
    Yb = matrixControlsNormalized,
    Z = Z,
    d.est = d.est,
    d.band = d.band,
    r = r,
    V.est.a = V.est, Tr.est.a = Tr.est,
    V.band.a = V.band, Tr.band.a = Tr.band,
    penalty = TRUE,
    lambda = lambda,
    alpha.grid = alpha.grid,
    adjust.sigma = TRUE
  )
  
  saveRDS(sccGroupComparison, file = sccComparisonPath)
  message("[INFO] SCC comparison result saved to ", sccComparisonPath)
  
} else {
  message("[INFO] SCC comparison result already exists — loading.")
  sccGroupComparison <- readRDS(sccComparisonPath)
}

# Custom plotting function

plotSCCcomparisonPanel <- function(scc,
                                   title = NULL,
                                   zlim = "auto",
                                   palette = "nih",
                                   label1 = "Group 1 Mean",
                                   label2 = "Group 2 Mean",
                                   label3 = "SCC Overlay") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required.")
  if (!requireNamespace("fields", quietly = TRUE)) stop("fields is required.")
  if (!requireNamespace("viridis", quietly = TRUE)) stop("viridis is required.")
  if (!requireNamespace("scales", quietly = TRUE)) stop("scales is required.")
  if (!requireNamespace("neuroSCC", quietly = TRUE)) stop("neuroSCC is required.")
  
  # 1. Extract coordinates and PET means
  Z <- scc$Z.band
  if (is.null(Z) || is.null(scc$Ya) || is.null(scc$Yb)) {
    stop("SCC object must include Z.band, Ya, and Yb.")
  }
  
  mean1_vals <- colMeans(scc$Ya)
  mean2_vals <- colMeans(scc$Yb)
  
  # Mask: only voxels with nonzero signal in at least one group
  mask <- (mean1_vals != 0 | mean2_vals != 0)
  
  df1 <- data.frame(x = Z[mask, 1], y = Z[mask, 2], value = mean1_vals[mask])
  df2 <- data.frame(x = Z[mask, 1], y = Z[mask, 2], value = mean2_vals[mask])
  df3 <- df2  # used for grayscale background
  
  # 2. Overlay points
  overlay <- neuroSCC::getPoints(scc)
  posPts <- overlay$positivePoints
  negPts <- overlay$negativePoints
  
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
  
  if (nrow(negPts) > 0) {
    p3 <- p3 + ggplot2::geom_point(data = negPts, ggplot2::aes(x = x, y = y), inherit.aes = FALSE,
                                   color = "red", size = 1.6, shape = 17)
  }
  if (nrow(posPts) > 0) {
    p3 <- p3 + ggplot2::geom_point(data = posPts, ggplot2::aes(x = x, y = y), inherit.aes = FALSE,
                                   color = "blue", size = 1.6, shape = 15)
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


# ---- Plot Group vs Group SCC Panel ----
plotSCCcomparisonPanel(
  scc = sccGroupComparison,
  label1 = "Alzheimer's Group Mean",
  label2 = "Control Group Mean",
  label3 = "Significant Regions",
  title  = "SCC Comparison: AD vs CN"
)

  
  
  
  
  
  