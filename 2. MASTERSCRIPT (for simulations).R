############################## ################### ############## ### 
##
## Script name: MASTERSCRIPT (for simulations)
##
## Purpose of script: Script for the calculation of sensibility, specificity, predictive value... 
## and other measures for comparing SCC and SPM using simulated data.
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
###                     1) PREAMBLE                      ###
### ==================================================== ###

#* Set working directory
setwd("~/GitHub/PhD-2023-SCC-vs-SPM-Group-vs-Group")

#* Options
options(scipen = 6, digits = 4)

#* Define axial slice to analyze
paramZ <- 35

#* Define CRAN packages
cran_packages <- c(
  "ggplot2", "patchwork", "fields", "viridis",     # Plotting
  "knitr", "kableExtra", "tibble", "magrittr",     # Tables & reports
  "tidyr", "dplyr", "scales",                      # Data wrangling
  "xtable", "neuroSCC"                             # Output and main package
)

#* Load CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

#* Load drat-hosted packages (non-CRAN)
drat_packages <- c("Triangulation", "ImageSCC", "BPST")
for (pkg in drat_packages) {
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

#* Clean up
rm(pkg, cran_packages, drat_packages)


### ==================================================== ###
###             2) CONTOURS OF NEURO-DATA               ###
### ==================================================== ###

#* Define file paths
mask_path <- "Auxiliary Files/new_mask.nii"
triangulation_path <- file.path("z35", "contour35.RData")

#* Only run if triangulation file does not exist
if (!file.exists(triangulation_path)) {
  
  message("[INFO] Generating triangulation from mask...")
  
  # Extract outer brain contour at level 0
  contours <- neuroSCC::neuroContour(mask_path, levels = 0)
  
  # Use external boundary only for triangulation
  triangulation <- Triangulation::TriMesh(contours[[1]], n = 15)
  
  # Save result for reuse
  if (!dir.exists("z35")) dir.create("z35", recursive = TRUE)
  save(triangulation, file = triangulation_path)
  
  message("[INFO] Triangulation saved to: ", triangulation_path)
  
} else {
  message("[INFO] Triangulation already exists. Skipping generation.")
}

### ==================================================== ###
###     3) CREATE SCC MATRIX FOR CONTROL GROUP          ###
### ==================================================== ###

#* Define output path
matrixCNPath <- file.path("z35", "SCC_CN.RData")

#* Only run if the matrix hasn't already been created
if (!file.exists(matrixCNPath)) {
  
  message("[INFO] Creating SCC matrix for Control group...")
  
  # Set working directory to folder containing input .nii.gz files
  setwd("PETimg_masked for simulations")
  
  # Create database for Controls (pattern: w00)
  databaseCN <- neuroSCC::databaseCreator(
    pattern = "_w00_",      # matches all control files
    control = TRUE,
    quiet = FALSE
  )
  
  # Return to root
  setwd("..")
  
  # Create subject-by-voxel matrix at z = paramZ
  matrixCN <- neuroSCC::matrixCreator(
    database = databaseCN,
    paramZ = paramZ,
    quiet = FALSE
  )
  
  # Mean normalization
  matrixCN <- neuroSCC::meanNormalization(matrixCN, quiet = FALSE)
  
  # Save for reuse
  save(matrixCN, file = matrixCNPath)
  message("[INFO] Control matrix saved to: ", matrixCNPath)
  
} else {
  message("[INFO] Control matrix already exists. Skipping.")
}

#* Clean up
rm(matrixCNPath, databaseCN, matrixCN)


### ==================================================== ###
### 4) CREATE SCC MATRICES FOR PATHOLOGICAL GROUP       ###
### ==================================================== ###

# Full region list in correct order (smallest to largest)
regions <- c("w32", "w79", "w214", "w271", "w413", "roiAD")
rois <- c(1, 2, 4, 6, 8)

# Set the folder where input NIfTI files are stored
inputDir <- "PETimg_masked for simulations"
outputDir <- "z35"

# Loop through region × ROI
for (region in regions) {
  for (roi in rois) {
    
    # Define expected output filename
    matrixADPath <- file.path(outputDir, paste0("SCC_", region, "_", roi, ".RData"))
    
    if (!file.exists(matrixADPath)) {
      
      message("[INFO] Processing region ", region, ", ROI ", roi, "...")
      
      # Change working directory to input folder
      setwd(inputDir)
      
      # Build regex pattern to match filenames
      pattern <- paste0("^masked_swwwC.*_", region, "_0_", roi, "_.*\\.nii$")
      
      # Create database for this region/ROI
      databaseAD <- neuroSCC::databaseCreator(
        pattern = pattern,
        control = FALSE,
        quiet = FALSE
      )
      
      # Return to root project folder
      setwd("..")
      
      # Create matrix at selected slice
      matrixAD <- neuroSCC::matrixCreator(
        databaseAD,
        paramZ = paramZ,
        quiet = FALSE
      )
      
      # Normalize the matrix
      matrixAD <- neuroSCC::meanNormalization(matrixAD, quiet = FALSE)
      
      # Save result
      save(matrixAD, file = matrixADPath)
      message("[INFO] Matrix saved to: ", matrixADPath)
      
    } else {
      message("[INFO] Matrix for region ", region, ", ROI ", roi, " already exists. Skipping.")
    }
  }
}

#* Clean up only path-related variables
rm(inputDir, outputDir, matrixADPath, pattern)


### ==================================================== ###
### 5) SCC ESTIMATION                                    ###
### ==================================================== ###

# Define inputs
triangulationPath <- file.path("z35", "contour35.RData")
matrixCNPath <- file.path("z35", "SCC_CN.RData")
resultsDir <- file.path("z35", "results")

# Load triangulation (VT)
load(triangulationPath)  # loads: VT
V.est <- as.matrix(VT[[1]])
Tr.est <- as.matrix(VT[[2]])
V.band <- V.est
Tr.band <- Tr.est

# Load control group matrix
load(matrixCNPath)  # loads: SCC_CN
matrixCN <- SCC_CN

# Define region/ROI grid
regions <- c("w32", "w79", "w214", "w271", "w413", "roiAD")
rois <- c(1, 2, 4, 6, 8)

# SCC hyperparameters
d.est <- 5
d.band <- 2
r <- 1
lambda <- 10^{seq(-6, 3, 0.5)}
alpha.grid <- c(0.10, 0.05, 0.01)

# Ensure results folder exists
if (!dir.exists(resultsDir)) dir.create(resultsDir, recursive = TRUE)

# Loop through all region × ROI
for (region in regions) {
  for (roi in rois) {
    
    resultPath <- file.path(resultsDir, paste0("SCC_COMP_", region, "_", roi, ".RData"))
    matrixADPath <- file.path("z35", paste0("SCC_", region, "_", roi, ".RData"))
    
    if (!file.exists(resultPath)) {
      
      message("[INFO] Running SCC estimation for ", region, " (ROI ", roi, ")...")
      
      load(matrixADPath)  # loads: matrixAD
      
      dims <- neuroSCC::getDimensions("Auxiliary Files/new_mask.nii")
      Z <- as.matrix(expand.grid(y = 1:dims$yDim, x = 1:dims$xDim)[, c(2, 1)])
      
      SCCcomp <- ImageSCC::scc.image(
        Ya = matrixAD,
        Yb = matrixCN,
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
      
      save(SCCcomp, file = resultPath)
      message("[INFO] SCC result saved to: ", resultPath)
      
      rm(matrixAD, SCCcomp)
      
    } else {
      message("[INFO] SCC result for ", region, " ROI ", roi, " already exists. Skipping.")
    }
  }
}

#* Clean up shared objects but keep matrixCN for next steps
rm(V.est, Tr.est, V.band, Tr.band, d.est, d.band, r, lambda, alpha.grid, resultPath, matrixADPath, triangulationPath, resultsDir)


### ==================================================== ###
### 6) SCC EVALUATION                                   ###
### ==================================================== ###

# Paths
resultsDir <- "z35/results"
roiDir <- "roisNormalizadas"
roiTablesDir <- file.path(roiDir, "tables")
maskPath <- "Auxiliary Files/new_mask.nii"

# Load coordinate grid from full image
dims <- neuroSCC::getDimensions(maskPath)

# Loop through region × ROI
for (region in regions) {
  for (roi in rois) {
    
    message("[INFO] Evaluating SCC result for region ", region, ", ROI ", roi, "...")
    
    # Fix naming discrepancy: use wroiAD only for ROI file input
    roiMaskRegion <- ifelse(region == "roiAD", "wroiAD", region)
    
    # Define SCC result path
    resultPath <- file.path(resultsDir, paste0("SCC_COMP_", region, "_", roi, ".RData"))
    outputDir <- file.path(resultsDir, paste0("ROI", roi))
    outputCSV <- file.path(outputDir, paste0("sens_esp_SCC_", region, "_", roi, ".csv"))
    
    # Ensure output folder exists
    if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
    
    # Load SCC result once
    load(resultPath)  # loads SCCcomp
    detectedPoints <- neuroSCC::getPoints(SCCcomp)$positivePoints
    
    # Temporary storage
    metricsList <- list()
    
    # Loop through subjects (C1 to C25)
    for (subject in 1:25) {
      
      subjectID <- paste0("C", subject)
      
      # Construct path to subject-specific ROI file
      roiFile <- file.path(
        roiDir,
        paste0("wwwx", roiMaskRegion, "_redim_crop_squ_flipLR_newDim_", subjectID, ".nii")
      )
      
      # Check ROI file exists
      if (!file.exists(roiFile)) {
        warning("[WARNING] ROI file not found: ", roiFile)
        next
      }
      
      # Process the ROI
      truePoints <- neuroSCC::processROIs(
        roiFile = roiFile,
        region = roiMaskRegion,
        number = subjectID,
        save = FALSE,
        verbose = FALSE
      )
      
      # Skip if no voxels at this slice
      sliceTruePoints <- subset(truePoints, z == paramZ & pet == 1)
      if (nrow(sliceTruePoints) == 0) {
        message("[WARNING] No ROI voxels at z = ", paramZ, " for ", subjectID, " — skipping.")
        next
      }
      
      # Calculate metrics
      subjectMetrics <- neuroSCC::calculateMetrics(
        detectedPoints = detectedPoints,
        truePoints = truePoints,
        totalCoords = dims,
        regionName = paste0(region, "_", roi)
      )
      
      subjectMetrics$subject <- subjectID
      metricsList[[length(metricsList) + 1]] <- subjectMetrics
    }
    
    # Combine all rows into one data frame
    if (length(metricsList) > 0) {
      allMetrics <- do.call(rbind, metricsList)
      readr::write_csv(allMetrics, outputCSV)
      message("[INFO] Metrics written to: ", outputCSV)
    } else {
      warning("[WARNING] No metrics computed for ", region, " ROI ", roi)
    }
    
    # Clean memory
    rm(SCCcomp, detectedPoints, metricsList, allMetrics)
  }
}

# Clean path-level objects
rm(resultsDir, roiDir, roiTablesDir, maskPath, outputDir, outputCSV, roiFile, resultPath)


### ==================================================== ###
### 6B) SPM EVALUATION                                  ###
### ==================================================== ###

# Paths
spmDir <- "z35/SPM"
roiDir <- "roisNormalizadas"
roiTablesDir <- file.path(roiDir, "tables")
resultsDir <- "z35/results"
maskPath <- "Auxiliary Files/new_mask.nii"

# Load coordinate grid
dims <- neuroSCC::getDimensions(maskPath)

# Loop through region × ROI
for (region in regions) {
  for (roi in rois) {
    
    message("[INFO] Evaluating SPM result for region ", region, ", ROI ", roi, "...")
    
    # Fix naming inconsistencies
    roiMaskRegion <- ifelse(region == "roiAD", "wroiAD", region)     # for ROI mask
    roiFolderRegion <- ifelse(region == "roiAD", "wroiAD", region)   # for SPM folder
    
    # Determine region index for correct folder lookup
    regionIndex <- which(regions == region)
    spmSubdir <- file.path(spmDir, paste0("ROI", regionIndex, "_", roiFolderRegion, "_0", roi))
    niftiFile <- file.path(spmSubdir, "binary.nii")
    
    # Result paths
    outputDir <- file.path(resultsDir, paste0("ROI", roi))
    outputCSV <- file.path(outputDir, paste0("sens_esp_SPM_", region, "_", roi, ".csv"))
    
    # Ensure output folder exists
    if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
    
    # Check for binary.nii
    if (!file.exists(niftiFile)) {
      warning("[WARNING] binary.nii not found at: ", niftiFile)
      next
    }
    
    # Load SPM-detected points (same for all subjects)
    detectedPoints <- neuroSCC::getSPMbinary(niftiFile, paramZ = paramZ)
    
    # Temporary storage for all subject evaluations
    metricsList <- list()
    
    for (subject in 1:25) {
      subjectID <- paste0("C", subject)
      
      # Construct subject-specific ROI file
      roiFile <- file.path(
        roiDir,
        paste0("wwwx", roiMaskRegion, "_redim_crop_squ_flipLR_newDim_", subjectID, ".nii")
      )
      
      # Check existence
      if (!file.exists(roiFile)) {
        warning("[WARNING] ROI file not found for subject ", subjectID, ": ", roiFile)
        next
      }
      
      # Process subject-specific ROI
      truePoints <- neuroSCC::processROIs(
        roiFile = roiFile,
        region = roiMaskRegion,
        number = subjectID,
        save = FALSE,
        verbose = FALSE
      )
      
      # Filter for z = paramZ and pet = 1
      sliceTruePoints <- subset(truePoints, z == paramZ & pet == 1)
      if (nrow(sliceTruePoints) == 0) {
        message("[WARNING] No ROI voxels at z = ", paramZ, " for ", subjectID, " — skipping.")
        next
      }
      
      # Evaluate
      subjectMetrics <- neuroSCC::calculateMetrics(
        detectedPoints = detectedPoints,
        truePoints = truePoints,
        totalCoords = dims,
        regionName = paste0(region, "_", roi)
      )
      
      subjectMetrics$subject <- subjectID
      metricsList[[length(metricsList) + 1]] <- subjectMetrics
    }
    
    # Write all evaluations to CSV
    if (length(metricsList) > 0) {
      allMetrics <- do.call(rbind, metricsList)
      readr::write_csv(allMetrics, outputCSV)
      message("[INFO] Metrics written to: ", outputCSV)
    } else {
      warning("[WARNING] No SPM metrics computed for ", region, " ROI ", roi)
    }
    
    # Clean
    rm(detectedPoints, metricsList, allMetrics)
  }
}

# Clean global temp vars
rm(spmDir, roiDir, roiTablesDir, resultsDir, maskPath, outputDir, outputCSV, roiFile, niftiFile, regionIndex)


### ==================================================== ###
### 7) SUMMARY TABLES                                   ###
### ==================================================== ###

library(dplyr)
library(readr)

# Storage for full dataset
SCC_vs_SPM_complete <- data.frame()

# Loop through region × ROI
for (region in regions) {
  for (roi in rois) {
    
    # Build folder path
    resultFolder <- file.path("z35", "results", paste0("ROI", roi))
    
    # Define file paths
    fileSCC <- file.path(resultFolder, paste0("sens_esp_SCC_", region, "_", roi, ".csv"))
    fileSPM <- file.path(resultFolder, paste0("sens_esp_SPM_", region, "_", roi, ".csv"))
    
    # If SCC exists, read and tag
    if (file.exists(fileSCC)) {
      tempSCC <- read_csv(fileSCC, show_col_types = FALSE)
      tempSCC <- tempSCC |>
        tidyr::separate(col = region, into = c("region", "roi"), sep = "_") |>
        dplyr::mutate(
          method = "SCC",
          roi = as.numeric(roi)
        ) |>
        dplyr::select(method, region, roi, dplyr::everything())
      SCC_vs_SPM_complete <- bind_rows(SCC_vs_SPM_complete, tempSCC)
    }
    
    # If SPM exists, read and tag
    if (file.exists(fileSPM)) {
      tempSPM <- read_csv(fileSPM, show_col_types = FALSE)
      tempSPM <- tempSPM |>
        tidyr::separate(col = region, into = c("region", "roi"), sep = "_") |>
        dplyr::mutate(
          method = "SPM",
          roi = as.numeric(roi)
        ) |>
        dplyr::select(method, region, roi, dplyr::everything())
      
      SCC_vs_SPM_complete <- dplyr::bind_rows(SCC_vs_SPM_complete, tempSPM)
    }
  }
}

# Create grouped summary (mean ± SD)
SCC_vs_SPM <- SCC_vs_SPM_complete |>
  group_by(method, region, roi) |>
  summarise(
    sensMEAN = mean(sensitivity, na.rm = TRUE),
    sensSD   = sd(sensitivity, na.rm = TRUE),
    espMEAN  = mean(specificity, na.rm = TRUE),
    espSD    = sd(specificity, na.rm = TRUE),
    ppvMEAN  = mean(PPV, na.rm = TRUE),
    ppvSD    = sd(PPV, na.rm = TRUE),
    npvMEAN  = mean(NPV, na.rm = TRUE),
    npvSD    = sd(NPV, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(method, region, roi)

# ====================================================
# === FILTER DATA FIRST (only selected regions & ROIs) ===
# ====================================================

# Keep only regions w32, w214, w271, roiAD
SCC_vs_SPM_complete <- SCC_vs_SPM_complete %>%
  filter(region %in% c("w32", "w214", "w271", "roiAD")) %>%
  filter(roi %in% c(1, 4, 8))  # keep hypo levels 1, 4, 8 only

SCC_vs_SPM <- SCC_vs_SPM %>%
  filter(region %in% c("w32", "w214", "w271", "roiAD")) %>%
  filter(roi %in% c(1, 4, 8))

# Update factor levels for new ROI labeling
SCC_vs_SPM_complete$region <- factor(SCC_vs_SPM_complete$region,
                       levels = c("w32", "w214", "w271", "roiAD"),
                       labels = c("ROI 1", "ROI 2", "ROI 3", "ROI 4"))
SCC_vs_SPM$region <- factor(SCC_vs_SPM$region,
                            levels = c("w32", "w214", "w271", "roiAD"),
                            labels = c("ROI 1", "ROI 2", "ROI 3", "ROI 4"))
SCC_vs_SPM_complete$roi <- factor(SCC_vs_SPM_complete$roi, levels = c(1, 4, 8), 
                                  labels = c("10", "40", "80"))
SCC_vs_SPM$roi <- factor(SCC_vs_SPM$roi, levels = c(1, 4, 8), labels = c("10", "40", "80"))

# Save summary version
saveRDS(SCC_vs_SPM, file = "z35/results/SCC_vs_SPM.RDS")
write_csv(SCC_vs_SPM, "z35/results/SCC_vs_SPM.csv")

# Save full version
saveRDS(SCC_vs_SPM_complete, file = "z35/results/SCC_vs_SPM_complete.RDS")
write_csv(SCC_vs_SPM_complete, "z35/results/SCC_vs_SPM_complete.csv")

# Clean up
rm(tempSCC, tempSPM, fileSCC, fileSPM, resultFolder)


### ==================================================== ###
### 7B) STATISTICAL SIGNIFICANCE TESTS (Paired T-Tests) ###
### ==================================================== ###

library(dplyr)
library(tidyr)
library(readr)

#* Helper to assign stars
get_stars <- function(p) {
  if (is.na(p)) return("")
  if (p <= 0.001) return("***")
  if (p <= 0.01)  return("**")
  if (p <= 0.05)  return("*")
  return("")
}

#* Ensure tibble
SCC_vs_SPM_complete <- as_tibble(SCC_vs_SPM_complete)

#* Initialize result list
pvalue_table <- list()

#* Loop over each region × roi
for (reg in levels(SCC_vs_SPM_complete$region)) {
  for (r in levels(SCC_vs_SPM_complete$roi)) {
    
    subset_data <- SCC_vs_SPM_complete %>%
      filter(region == reg, roi == r)
    
    # Check both methods are present
    if (!all(c("SCC", "SPM") %in% unique(subset_data$method))) next
    
    # Prepare data for paired comparison
    wide_data <- subset_data %>%
      dplyr::select(method, subject, sensitivity, specificity, PPV, NPV) %>%
      pivot_wider(names_from = method, values_from = c(sensitivity, specificity, PPV, NPV)) %>%
      drop_na()
    
    # Run paired t-tests
    t_sens <- t.test(wide_data$sensitivity_SCC, wide_data$sensitivity_SPM, paired = TRUE)
    t_esp  <- t.test(wide_data$specificity_SCC, wide_data$specificity_SPM, paired = TRUE)
    t_ppv  <- t.test(wide_data$PPV_SCC, wide_data$PPV_SPM, paired = TRUE)
    t_npv  <- t.test(wide_data$NPV_SCC, wide_data$NPV_SPM, paired = TRUE)
    
    # Store results
    pvalue_table[[paste(reg, r, sep = "_")]] <- tibble(
      region = reg,
      roi = r,
      p_sens = t_sens$p.value,
      sig_sens = get_stars(t_sens$p.value),
      p_esp = t_esp$p.value,
      sig_esp = get_stars(t_esp$p.value),
      p_ppv = t_ppv$p.value,
      sig_ppv = get_stars(t_ppv$p.value),
      p_npv = t_npv$p.value,
      sig_npv = get_stars(t_npv$p.value)
    )
  }
}

#* Combine into one table
pvalue_table_Groups <- bind_rows(pvalue_table)

#* Save results
write_csv(pvalue_table_Groups, "z35/results/pvalue_table_group.csv")
saveRDS(pvalue_table_Groups, "z35/results/pvalue_table_group.RDS")

#* Print preview
print(pvalue_table_Groups)

#* Clean up
rm(get_stars, wide_data, subset_data, t_sens, t_esp, t_ppv, t_npv)


### ==================================================== ###
### 8) VISUALIZATIONS                                   ###
### ==================================================== ###

# Load data (summary and full evaluation)
referencia <- readRDS(paste0("~/GitHub/PhD-2023-SCC-vs-SPM-Group-vs-Group/z", 
                             as.numeric(paramZ), "/results/SCC_vs_SPM.RDS"))
table <- readRDS(paste0("~/GitHub/PhD-2023-SCC-vs-SPM-Group-vs-Group/z", 
                        as.numeric(paramZ), "/results/SCC_vs_SPM_complete.RDS"))

# Load necessary packages
library(tidyverse)
library(lemon)
library(gridExtra)
library(ggridges)
library(viridis)
library(dotwhisker)

# Set save directory
setwd(paste0("~/GitHub/PhD-2023-SCC-vs-SPM-Group-vs-Group/z", as.numeric(paramZ), "/Figures"))

### ----------------------------------------------------------------
### Sensitivity, Specificity, PPV, NPV across all filtered regions
### ----------------------------------------------------------------

plot_metric_by_region <- function(metric_name, y_label) {
  ggplot(data = table, aes(x = roi, y = .data[[metric_name]])) +
    geom_boxplot(aes(fill = method), outlier.size = 0.5, lwd = 0.25) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0, 100)) +
    xlab("Hypoactivity (%)") +
    ylab(y_label) +
    facet_wrap(~region, ncol = 2) +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal(base_family = "serif") +
    theme(panel.border = element_blank(),
          axis.line = element_line(),
          legend.text = element_text(size = 14))
}

graph_sens <- plot_metric_by_region("sensitivity", "Sensitivity (%)")
graph_esp  <- plot_metric_by_region("specificity", "Specificity (%)")
graph_ppv  <- plot_metric_by_region("PPV", "Positive Predictive Value (%)")
graph_npv  <- plot_metric_by_region("NPV", "Negative Predictive Value (%)")

ggsave("sens_FILTERED.png", plot = graph_sens, width = 28.95, height = 18.3, units = "cm", dpi = 600)
ggsave("esp_FILTERED.png", plot = graph_esp, width = 28.95, height = 18.3, units = "cm", dpi = 600)
ggsave("ppv_FILTERED.png", plot = graph_ppv, width = 28.95, height = 18.3, units = "cm", dpi = 600)
ggsave("npv_FILTERED.png", plot = graph_npv, width = 28.95, height = 18.3, units = "cm", dpi = 600)

combined_filtered <- grid.arrange(graph_sens, graph_esp, graph_ppv, graph_npv, ncol = 2)
ggsave("summary_metrics_FILTERED.png", plot = combined_filtered, width = 38, height = 28, units = "cm", dpi = 600)

### ==================================================== ###
### OPTIONAL VISUALIZATION EXPERIMENTS (Ridge + Heatmap) ###
### ==================================================== ###

### -----------------------------------------------
### Ridge plot: Sensitivity by Region
### -----------------------------------------------

ridge_plot_sens <- ggplot(table, aes(x = sensitivity, y = region, fill = method)) +
  geom_density_ridges(alpha = 0.6, scale = 1.2, rel_min_height = 0.01) +
  scale_fill_brewer(palette = "Set1") +
  theme_ridges(font_family = "serif") +
  labs(x = "Sensitivity (%)", y = "Region", title = "Density of Sensitivity by Region")

# ggsave("ridge_sensitivity_FILTERED.png", ridge_plot_sens, width = 28, height = 18, units = "cm", dpi = 600)

### -----------------------------------------------
### Ridge plot: PPV by Region
### -----------------------------------------------

ridge_plot_ppv <- ggplot(table, aes(x = PPV, y = region, fill = method)) +
  geom_density_ridges(alpha = 0.6, scale = 1.2, rel_min_height = 0.01) +
  scale_fill_brewer(palette = "Set1") +
  theme_ridges(font_family = "serif") +
  labs(x = "Positive Predictive Value (%)", y = "Region", title = "Density of PPV by Region")

# ggsave("ridge_ppv_FILTERED.png", ridge_plot_ppv, width = 28, height = 18, units = "cm", dpi = 600)

### -----------------------------------------------
### Double-faceted Heatmap: Sensitivity
### -----------------------------------------------

heatmap_sens_facet <- ggplot(referencia, aes(x = factor(roi), y = region, fill = sensMEAN)) +
  geom_tile(color = "white") +
  facet_wrap(~method) +
  scale_fill_viridis(name = "Sensitivity", option = "C", limits = c(0, 100)) +
  labs(x = "Hypoactivity Level (%)", y = "Region", title = "Sensitivity (mean) by Method") +
  theme_minimal(base_family = "serif")

# ggsave("heatmap_sensitivity_FILTERED.png", heatmap_sens_facet, width = 28, height = 14, units = "cm", dpi = 600)

### -----------------------------------------------
### Double-faceted Heatmap: PPV
### -----------------------------------------------

heatmap_ppv_facet <- ggplot(referencia, aes(x = factor(roi), y = region, fill = ppvMEAN)) +
  geom_tile(color = "white") +
  facet_wrap(~method) +
  scale_fill_viridis(name = "PPV", option = "C", limits = c(0, 100)) +
  labs(x = "Hypoactivity Level (%)", y = "Region", title = "PPV (mean) by Method") +
  theme_minimal(base_family = "serif")

# ggsave("heatmap_ppv_FILTERED.png", heatmap_ppv_facet, width = 28, height = 14, units = "cm", dpi = 600)


### ==================================================== ###
### 9) BRAIN-LEVEL VISUALIZATIONS                        ###
### ==================================================== ###

# Load custom visualization helpers
setwd("~/GitHub/PhD-2023-SCC-vs-SPM-Group-vs-Group")
source("Visualization Functions.R")

template <- system.file("extdata", "syntheticControl1.nii.gz", package = "neuroSCC")
dims <- neuroSCC::getDimensions(template)
Z <- as.matrix(expand.grid(y = 1:dims$yDim, x = 1:dims$xDim)[, c(2, 1)])

# Load SCC objects if not already loaded
load("~/GitHub/PhD-2023-SCC-vs-SPM-Group-vs-Group/z35/SCC_CN.RData")
load("z35/SCC_roiAD_1.RData")
load("z35/SCC_w32_4.RData")
load("z35/results/SCC_COMP_w32_1.RData")
load("z35/results/SCC_COMP_w214_4.RData") # example alternative
sccOneGroup <- readRDS("~/GitHub/PhD-2023-SCC-vs-SPM-Group-vs-Group/z35/sccOneGroup.RDS")

# Set working directory for saving figures
setwd(paste0("~/GitHub/PhD-2023-SCC-vs-SPM-Group-vs-Group/z", as.numeric(paramZ), "/Figures"))

# Example 1: Visualize SCC for a Control group (one-group SCC)

oneGroupSCC <- readRDS("~/GitHub/PhD-2023-SCC-vs-SPM-Group-vs-Group/oneGroupSCC.rds")

plotSCCpanel(
  scc = oneGroupSCC,                     # Alternative: SCC_AD
  title = "",   # Try: "SCC FOR PATHOLOGICAL GROUP"
  zlim = "auto",                     # Try: c(-0.5, 0.5)
  palette = "nih"                    # Try: "viridis", "gray"
)

# Example 2: Visualize Group vs Group SCC Comparison
plotSCCcomparisonPanel(
  scc = SCC_COMP,              # Alternative: SCC_COMP_w214_4
  title = "CONTROL vs PATHOLOGICAL COMPARISON",
  label1 = "Pathological Group Mean",   # Can swap label order depending on Ya, Yb
  label2 = "Control Group Mean",
  label3 = "Detected Significant Regions",
  zlim = "auto",                        # Try: c(-1, 1)
  palette = "nih"                       # Try: "viridis"
)

# Alternative Group Comparison (example commented)
# plotSCCcomparisonPanel(
#   scc = SCC_COMP_w214_4,
#   title = "Another Group Comparison",
#   label1 = "AD Group Mean",
#   label2 = "Control Group Mean",
#   label3 = "SCC-Detected Hypo/Hyperactivity"
# )

# Example 3: Validation Panel (ROI vs SCC vs SPM points)
# Precompute supporting objects:

# Load the background template and control matrix if needed
roiFile <- "Auxiliary Files/new_mask.nii"  # or another mask if you want
load("~/GitHub/PhD-2023-SCC-vs-SPM-Group-vs-Group/z35/SCC_CN.RData")    # or load directly if not present

# Load ROI points
roiPoints <- neuroSCC::processROIs(
  roiFile = file.path("roisNormalizadas", "wwwxw32_redim_crop_squ_flipLR_newDim_C1.nii"),
  region = "w32",
  number = "1",
  save = FALSE,
  verbose = FALSE
)
roiPoints <- subset(roiPoints, z == paramZ & pet == 1, select = c("x", "y"))

# Load SCC comparison
load("z35/results/SCC_COMP_w32_1.RData")
sccPoints <- neuroSCC::getPoints(SCC_COMP)

# Load SPM detected points
spmBinary <- "z35/SPM/ROI1_w32_01/binary.nii"
spmPoints <- neuroSCC::getSPMbinary(spmBinary, paramZ = paramZ)

# Plot
plotValidationPanel(
  template = roiFile,
  backgroundMatrix = SCC_CN,
  roiPoints = roiPoints,
  sccPoints = sccPoints,
  spmPoints = spmPoints,
  title = "Performance Validation Panel",
  label1 = "Ground Truth (ROI)",
  label2 = "SCC Detected",
  label3 = "SPM Detected"
)

# Alternative Validation:
# Use SCC_COMP_w214_4 and corresponding ROI + SPM binary

















