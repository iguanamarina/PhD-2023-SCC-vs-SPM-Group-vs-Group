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
### 6A) SCC EVALUATION                                   ###
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
### 7B) STATISTICAL SIGNIFICANCE TESTS with testCompareR ###
### ==================================================== ###

# NOTE for future users:
# In this section we pivot to a more robust and appropriate method for
# estimating sensitivity, specificity, PPV, NPV, and likelihood ratios — and 
# comparing SCC vs. SPM using the `compareR()` function from the testCompareR package.
#
# This replaces older paired t-tests which are not appropriate for proportion-based metrics.
# For theoretical background, see Roldán-Nofuentes (2020), BMC Med Res Methodol.
# DOI: 10.1186/s12874-020-00988-y

# Load required packages
if (!requireNamespace("testCompareR", quietly = TRUE)) {
  install.packages("testCompareR")
}
library(testCompareR)
library(readr)
library(dplyr)

# Load custom triplet builder
source("~/GitHub/PhD-2023-SCC-vs-SPM-Group-vs-Group/Contrastes de Hipótesis/generateCompareRTriplets.R")

# Define key paths and voxel grid
maskPath <- "Auxiliary Files/new_mask.nii"
dims <- neuroSCC::getDimensions(maskPath)
grid <- expand.grid(y = 1:dims$yDim, x = 1:dims$xDim)[, c("x", "y")]

# Full unfiltered region and ROI list
regions <- c("w32", "w79", "w214", "w271", "w413", "roiAD")
rois <- c(1, 2, 4, 6, 8)

# Initialize result storage
pvalResults <- list()

# Region × ROI loop
for (region in regions) {
  for (roi in rois) {
    
    message("[INFO] Running compareR() for region ", region, ", ROI ", roi, "...")
    
    # Adjust naming for ROI masks
    roiMaskRegion <- ifelse(region == "roiAD", "wroiAD", region)
    
    # Load SCC results
    sccPath <- file.path("z35/results", paste0("SCC_COMP_", region, "_", roi, ".RData"))
    if (!file.exists(sccPath)) {
      warning("[WARNING] Missing SCC_COMP file: ", sccPath)
      next
    }
    load(sccPath)  # loads SCC_COMP
    sccPoints <- neuroSCC::getPoints(SCC_COMP)$positivePoints
    
    # Load SPM results
    regionIndex <- match(region, c("w32", "w79", "w214", "w271", "w413", "roiAD"))
    spmDir <- file.path("z35/SPM", paste0("ROI", regionIndex, "_", roiMaskRegion, "_0", roi))
    spmFile <- file.path(spmDir, "binary.nii")
    if (!file.exists(spmFile)) {
      warning("[WARNING] Missing SPM binary.nii file: ", spmFile)
      next
    }
    spmPoints <- neuroSCC::getSPMbinary(spmFile, paramZ = paramZ)
    
    # Collect triplets across all subjects
    tripletList <- list()
    for (subject in 1:25) {
      subjectID <- paste0("C", subject)
      roiFile <- file.path(
        "roisNormalizadas",
        paste0("wwwx", roiMaskRegion, "_redim_crop_squ_flipLR_newDim_", subjectID, ".nii")
      )
      if (!file.exists(roiFile)) {
        warning("[WARNING] ROI file not found: ", roiFile)
        next
      }
      
      truePoints <- neuroSCC::processROIs(
        roiFile = roiFile,
        region = roiMaskRegion,
        number = subjectID,
        save = FALSE,
        verbose = FALSE
      )
      sliceTruePoints <- subset(truePoints, z == paramZ & pet == 1, select = c("x", "y"))
      if (nrow(sliceTruePoints) == 0) {
        message("[WARNING] No ROI voxels at z = ", paramZ, " for ", subjectID, " — skipping.")
        next
      }
      
      triplet <- generateCompareRTriplets(
        grid = grid,
        sccCoords = sccPoints,
        spmCoords = spmPoints,
        roiCoords = sliceTruePoints
      )
      tripletList[[length(tripletList) + 1]] <- triplet
    }
    
    # Merge triplets and convert to numeric
    if (length(tripletList) == 0) {
      warning("[WARNING] No valid subjects for ", region, " ROI ", roi)
      next
    }
    
    allTriplets <- bind_rows(tripletList) %>%
      mutate(across(everything(), as.numeric))
    
    # Run compareR
    result <- compareR(
      df = allTriplets,
      test1 = "SCC_pred",
      test2 = "SPM_pred",
      gold = "ROI_truth",
      interpret = FALSE,
      multi_corr = "holm",
      alpha = 0.05,
      sesp = TRUE,
      ppvnpv = TRUE,
      plrnlr = TRUE,  # enable LR+ and LR−
      test.names = c("SCC", "SPM"),
      dp = 2 # number of decimals
    )
    
    # Assemble row with full output
    row <- tibble(
      region = region,
      roi = roi,
      
      # Sensitivity
      sens_SCC = result$acc$accuracies$SCC["Sensitivity", "Estimate"],
      se_sens_SCC = result$acc$accuracies$SCC["Sensitivity", "SE"],
      sens_SPM = result$acc$accuracies$SPM["Sensitivity", "Estimate"],
      se_sens_SPM = result$acc$accuracies$SPM["Sensitivity", "SE"],
      p_sens = result$acc$sens.p.adj,
      
      # Specificity
      spec_SCC = result$acc$accuracies$SCC["Specificity", "Estimate"],
      se_spec_SCC = result$acc$accuracies$SCC["Specificity", "SE"],
      spec_SPM = result$acc$accuracies$SPM["Specificity", "Estimate"],
      se_spec_SPM = result$acc$accuracies$SPM["Specificity", "SE"],
      p_spec = result$acc$spec.p.adj,
      
      # PPV
      ppv_SCC = result$pv$predictive.values$SCC["PPV", "Estimate"],
      se_ppv_SCC = result$pv$predictive.values$SCC["PPV", "SE"],
      ppv_SPM = result$pv$predictive.values$SPM["PPV", "Estimate"],
      se_ppv_SPM = result$pv$predictive.values$SPM["PPV", "SE"],
      p_ppv = result$pv$ppv.p.adj,
      
      # NPV
      npv_SCC = result$pv$predictive.values$SCC["NPV", "Estimate"],
      se_npv_SCC = result$pv$predictive.values$SCC["NPV", "SE"],
      npv_SPM = result$pv$predictive.values$SPM["NPV", "Estimate"],
      se_npv_SPM = result$pv$predictive.values$SPM["NPV", "SE"],
      p_npv = result$pv$npv.p.adj,
      
      # Likelihood Ratios (PLR = LR+, NLR = LR−)
      lrpos_SCC = result$lr$likelihood.ratios$SCC["PLR", "Estimate"],
      se_lrpos_SCC = result$lr$likelihood.ratios$SCC["PLR", "SE"],
      lrpos_SPM = result$lr$likelihood.ratios$SPM["PLR", "Estimate"],
      se_lrpos_SPM = result$lr$likelihood.ratios$SPM["PLR", "SE"],
      p_lrpos = result$lr$plr.p.adj,
      
      lrneg_SCC = result$lr$likelihood.ratios$SCC["NLR", "Estimate"],
      se_lrneg_SCC = result$lr$likelihood.ratios$SCC["NLR", "SE"],
      lrneg_SPM = result$lr$likelihood.ratios$SPM["NLR", "Estimate"],
      se_lrneg_SPM = result$lr$likelihood.ratios$SPM["NLR", "SE"],
      p_lrneg = result$lr$nlr.p.adj
      
    )
    
    # Store result
    pvalResults[[length(pvalResults) + 1]] <- row
  }
}

# Combine into table
pvalTable <- bind_rows(pvalResults)

# Filter to final subset
pvalTable <- pvalTable %>%
  filter(region %in% c("w32", "w214", "w271", "roiAD"),
         roi %in% c(1, 4, 8))

# Significance stars
getStars <- function(p) {
  if (is.na(p)) return("")
  if (p <= 0.001) return("***")
  if (p <= 0.01)  return("**")
  if (p <= 0.05)  return("*")
  return("")
}

# Apply significance stars
pvalTable <- pvalTable %>%
  mutate(
    sig_sens   = sapply(p_sens, getStars),
    sig_spec   = sapply(p_spec, getStars),
    sig_ppv    = sapply(p_ppv, getStars),
    sig_npv    = sapply(p_npv, getStars),
    sig_lrpos  = sapply(p_lrpos, getStars),
    sig_lrneg  = sapply(p_lrneg, getStars)
  )

# Reorder by metric
pvalTable <- pvalTable %>%
  dplyr::select(
    region, roi,
    
    # Sensibilidad
    sens_SCC, se_sens_SCC, sens_SPM, se_sens_SPM, p_sens, sig_sens,
    
    # Especificidad
    spec_SCC, se_spec_SCC, spec_SPM, se_spec_SPM, p_spec, sig_spec,
    
    # Valor Predictivo Positivo
    ppv_SCC, se_ppv_SCC, ppv_SPM, se_ppv_SPM, p_ppv, sig_ppv,
    
    # Valor Predictivo Negativo
    npv_SCC, se_npv_SCC, npv_SPM, se_npv_SPM, p_npv, sig_npv,
    
    # Likelihood Ratio Positivo
    lrpos_SCC, se_lrpos_SCC, lrpos_SPM, se_lrpos_SPM, p_lrpos, sig_lrpos,
    
    # Likelihood Ratio Negativo
    lrneg_SCC, se_lrneg_SCC, lrneg_SPM, se_lrneg_SPM, p_lrneg, sig_lrneg
  )

# Save results
write_csv(pvalTable, "z35/results/pvalue_table_compareR.csv")
saveRDS(pvalTable, "z35/results/pvalue_table_compareR.RDS")
pvalue_table_compareR <- readRDS("~/GitHub/PhD-2023-SCC-vs-SPM-Group-vs-Group/z35/results/pvalue_table_compareR.RDS")

# Función para imprimir resumen de resultados para una fila de pvalTable
print_summary <- function(data, region_value, roi_value) {
  row <- data %>%
    filter(region == region_value, roi == roi_value) %>%
    slice(1)
  
  cat(glue::glue("
Región {region_value}, ROI {roi_value}:

- Sensibilidad:
  SCC = {round(row$sens_SCC, 2)}% ± {round(row$se_sens_SCC, 2)} {ifelse(row$sig_sens != '', paste0('(', row$sig_sens, ')'), '(sin significancia)')}
  SPM = {round(row$sens_SPM, 2)}% ± {round(row$se_sens_SPM, 2)}

- Especificidad:
  SCC = {round(row$spec_SCC, 2)}% ± {round(row$se_spec_SCC, 2)} (sin p-valor)
  SPM = {round(row$spec_SPM, 2)}% ± {round(row$se_spec_SPM, 2)}

- PPV:
  SCC = {round(row$ppv_SCC, 2)}% ± {round(row$se_ppv_SCC, 2)} {ifelse(row$sig_ppv != '', paste0('(', row$sig_ppv, ')'), '(sin significancia)')}
  SPM = {round(row$ppv_SPM, 2)}% ± {round(row$se_ppv_SPM, 2)}

- NPV:
  SCC = {round(row$npv_SCC, 2)}% ± {round(row$se_npv_SCC, 2)} {ifelse(row$sig_npv != '', paste0('(', row$sig_npv, ')'), '(sin significancia)')}
  SPM = {round(row$npv_SPM, 2)}% ± {round(row$se_npv_SPM, 2)}

"))
}

# Imprimir todas las combinaciones
for (i in seq_len(nrow(pvalue_table_compareR))) {
  print_summary(
    data = pvalue_table_compareR,
    region_value = pvalue_table_compareR$region[i],
    roi_value = pvalue_table_compareR$roi[i]
  )
}


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
### Sensitivity, Specificity, PPV, NPV across all filtered regions with asterisks
### ----------------------------------------------------------------

plot_roi_grouped_with_asterisks <- function(table, pval_table, p_col = "p_sens", metric = "sensitivity", y_label = "Sensitivity (%)") {
  library(ggplot2)
  library(dplyr)
  
  dodge_width <- 0.8
  
  # 1. Prepare the input table
  table <- table %>%
    mutate(
      method = factor(method, levels = c("SCC", "SPM")),
      roi = factor(roi, levels = c("10", "40", "80")),
      region = factor(region)
    )
  
  # 2. Clean and label brackets (only asterisks)
  pval_clean <- pval_table %>%
    dplyr::select(region, roi, !!sym(p_col)) %>%
    mutate(
      region = as.character(region),
      roi = as.character(roi),
      group1 = "SCC",
      group2 = "SPM",
      label = case_when(
        !!sym(p_col) <= 0.001 ~ "***",
        !!sym(p_col) <= 0.01  ~ "**",
        !!sym(p_col) <= 0.05  ~ "*",
        TRUE                  ~ "ns"
      )
    )
  
  # 3. Match panel structure
  valid_panels <- table %>%
    distinct(region, roi) %>%
    mutate(across(everything(), as.character))
  pval_clean <- semi_join(pval_clean, valid_panels, by = c("region", "roi"))
  
  # 4. Compute height of bracket
  bracket_y <- table %>%
    group_by(region, roi) %>%
    summarise(y.position = max(.data[[metric]], na.rm = TRUE) + 5, .groups = "drop") %>%
    mutate(across(c(region, roi), as.character))
  
  # 5. Merge brackets
  bracket_data <- pval_clean %>%
    left_join(bracket_y, by = c("region", "roi")) %>%
    mutate(
      x = as.numeric(factor(roi, levels = c("10", "40", "80"))),
      x1 = x - dodge_width / 4,
      x2 = x + dodge_width / 4,
      label_y = y.position + 2
    )
  
  # 6. Final plot
  ggplot(table, aes(x = roi, y = .data[[metric]], fill = method)) +
    geom_boxplot(
      position = position_dodge(width = dodge_width),
      width = 0.6,
      outlier.size = 0.5,
      lwd = 0.25
    ) +
    facet_wrap(~region, ncol = 2) +
    scale_fill_brewer(palette = "Set1") +
    scale_y_continuous(
      limits = c(0, NA),
      breaks = c(25, 50, 75, 100),
      expand = expansion(mult = c(0, 0.1))
    ) +
    xlab("Hypoactivity (%)") + ylab(y_label) +
    theme_minimal(base_family = "serif") +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(),
      legend.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 15),
      panel.spacing = unit(1, "lines")
    ) +
    # Bracket bar
    geom_segment(
      data = bracket_data,
      aes(x = x1, xend = x2, y = y.position, yend = y.position),
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    # Left vertical arm
    geom_segment(
      data = bracket_data,
      aes(x = x1, xend = x1, y = y.position, yend = y.position - 2),
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    # Right vertical arm
    geom_segment(
      data = bracket_data,
      aes(x = x2, xend = x2, y = y.position, yend = y.position - 2),
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    # Asterisk-only label
    geom_text(
      data = bracket_data,
      aes(x = x, y = label_y, label = label),
      size = 5,
      family = "serif",
      inherit.aes = FALSE
    )
}

# Set wd() for figure export
setwd("~/GitHub/PhD-2023-SCC-vs-SPM-Group-vs-Group/z35/Figures")

# Create all 4 plots
graph_sens <- plot_roi_grouped_with_asterisks(
  table = table,
  pval_table = pvalue_table_Groups,
  p_col = "p_sens",
  metric = "sensitivity",
  y_label = "Sensitivity (%)"
)

graph_esp <- plot_roi_grouped_with_asterisks(
  table = table,
  pval_table = pvalue_table_Groups,
  p_col = "p_esp",
  metric = "specificity",
  y_label = "Specificity (%)"
)

graph_ppv <- plot_roi_grouped_with_asterisks(
  table = table,
  pval_table = pvalue_table_Groups,
  p_col = "p_ppv",
  metric = "PPV",
  y_label = "Positive Predictive Value (%)"
)

graph_npv <- plot_roi_grouped_with_asterisks(
  table = table,
  pval_table = pvalue_table_Groups,
  p_col = "p_npv",
  metric = "NPV",
  y_label = "Negative Predictive Value (%)"
)

# Save to PNGs
png("sens_FILTERED.png", width = 2895, height = 1830, res = 300)
print(graph_sens)
dev.off()

# Specificity
png("esp_FILTERED.png", width = 2895, height = 1830, res = 300)
print(graph_esp)
dev.off()

# PPV
png("ppv_FILTERED.png", width = 2895, height = 1830, res = 300)
print(graph_ppv)
dev.off()

# NPV
png("npv_FILTERED.png", width = 2895, height = 1830, res = 300)
print(graph_npv)
dev.off()


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
### Double-faceted Heatmap: Sensitivity & Specificity
### -----------------------------------------------

library(ggplot2)
library(patchwork)
library(grid)  # for unit()

# Heatmap for Sensitivity
heatmap_sens_facet <- ggplot(referencia, aes(x = factor(roi), y = region, fill = sensMEAN)) +
  geom_tile(color = "white") +
  facet_wrap(~method) +
  scale_fill_gradient2(
    name = NULL,
    low = "red", mid = "yellow", high = "green",
    midpoint = 50,
    limits = c(0, 100)
  ) +
  labs(
    x = "Hypoactivity Level (%)",
    y = "Region",
    title = "Mean Sensitivity by Method"
  ) +
  theme_minimal(base_family = "serif") +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(2, "cm")
  )

# Heatmap for Specificity
heatmap_spec_facet <- ggplot(referencia, aes(x = factor(roi), y = region, fill = espMEAN)) +
  geom_tile(color = "white") +
  facet_wrap(~method) +
  scale_fill_gradient2(
    name = NULL,
    low = "red", mid = "yellow", high = "green",
    midpoint = 50,
    limits = c(0, 100)
  ) +
  labs(
    x = "Hypoactivity Level (%)",
    y = "",  # preserve vertical space
    title = "Mean Specificity by Method"
  ) +
  theme_minimal(base_family = "serif") +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(2, "cm")
  )

# Combine
combined_heatmap_sens_esp <- heatmap_sens_facet + heatmap_spec_facet +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )

# Show
combined_heatmap_sens_esp

# Save
ggsave("combined_heatmap_sens_esp.png", combined_heatmap_sens_esp, width = 28, height = 20, units = "cm", dpi = 600)


### ------------------------------------------------------
### Double-faceted Heatmap: PPV & NPV
### ------------------------------------------------------

# Heatmap for PPV
heatmap_ppv_facet <- ggplot(referencia, aes(x = factor(roi), y = region, fill = ppvMEAN)) +
  geom_tile(color = "white") +
  facet_wrap(~method) +
  scale_fill_gradient2(
    name = NULL,
    low = "red", mid = "yellow", high = "green",
    midpoint = 50,
    limits = c(0, 100)
  ) +
  labs(
    x = "Hypoactivity Level (%)",
    y = "Region",
    title = "Mean PPV by Method"
  ) +
  theme_minimal(base_family = "serif") +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(2, "cm")
  )

# Heatmap for NPV
heatmap_npv_facet <- ggplot(referencia, aes(x = factor(roi), y = region, fill = npvMEAN)) +
  geom_tile(color = "white") +
  facet_wrap(~method) +
  scale_fill_gradient2(
    name = NULL,
    low = "red", mid = "yellow", high = "green",
    midpoint = 50,
    limits = c(0, 100)
  ) +
  labs(
    x = "Hypoactivity Level (%)",
    y = "",  # preserve space, blank label
    title = "Mean NPV by Method"
  ) +
  theme_minimal(base_family = "serif") +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(2, "cm")
  )

# Combine
combined_heatmap_ppv_npv <- heatmap_ppv_facet + heatmap_npv_facet +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )

# Show
combined_heatmap_ppv_npv

# Save
ggsave("combined_heatmap_ppv_npv.png", combined_heatmap_ppv_npv, width = 28, height = 20, units = "cm", dpi = 600)



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
  roiFile = file.path("roisNormalizadas", "wwwxwroiAD_redim_crop_squ_flipLR_newDim_C2.nii"),
  region = "roiAD",
  number = "1",
  save = FALSE,
  verbose = FALSE
)
roiPoints <- subset(roiPoints, z == paramZ & pet == 1, select = c("x", "y"))

# Load SCC comparison
load("z35/results/SCC_COMP_roiAD_8.RData")
sccPoints <- neuroSCC::getPoints(SCC_COMP)

# Load SPM detected points
spmBinary <- "z35/SPM/ROI6_wroiAD_08/binary.nii"
spmPoints <- neuroSCC::getSPMbinary(spmBinary, paramZ = paramZ)

# Plot
plotValidationPanel(
  template = roiFile,
  backgroundMatrix = SCC_CN,
  roiPoints = roiPoints,
  sccPoints = sccPoints,
  spmPoints = spmPoints,
  title = "",
  label1 = "Ground Truth (ROI)",
  label2 = "SCC Detected",
  label3 = "SPM Detected"
)

# Alternative Validation:
# Use SCC_COMP_w214_4 and corresponding ROI + SPM binary

















