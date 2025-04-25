### ================================================================ ###
###     SCC STABILITY AUDIT ACROSS HYPOACTIVITY LEVELS (ROI %)       ###
### ================================================================ ###
# This script verifies that different ROI-level simulations (e.g., 10% vs 80%)
# are being properly used in SCC, and explains why the results may still be identical.
# It assumes files are located in z35/ and z35/results.

library(dplyr)
library(ggplot2)

# Setup parameters
region <- "w32"              # region to audit
roiLevels <- c(1, 2, 4, 6, 8) # levels of induced hypoactivity
subjectID <- "C10"           # a subject included in the analysis
paramZ <- 35
xDim <- 91

# 1. Load all SCC result files and their corresponding matrixAD inputs

basePath <- "z35/results"
matrixADs <- list()
sccObjects <- list()
meanCurves <- list()

for (roi in roiLevels) {
  roi_chr <- as.character(roi)
  
  # Load SCC result
  sccPath <- file.path(basePath, paste0("SCC_COMP_", region, "_", roi, ".RData"))
  load(sccPath)  # loads SCCcomp
  sccObjects[[roi_chr]] <- SCCcomp
  meanCurves[[roi_chr]] <- SCCcomp$scc[, 1, 3]  # estimated mean function
  
  # Load corresponding matrixAD used in SCC
  matrixPath <- file.path("z35", paste0("SCC_", region, "_", roi, ".RData"))
  tempEnv <- new.env()
  load(matrixPath, envir = tempEnv)
  for (objName in ls(tempEnv)) {
    obj <- tempEnv[[objName]]
    if (is.matrix(obj)) {
      matrixADs[[roi_chr]] <- obj
      cat("✔️  Loaded matrix for ROI", roi, " → object:", objName, "\n")
      break
    }
  }
}

# 2. Check if any matrix inputs were reused by mistake

cat("\n❓ Are any input matrices (matrixADs) identical?\n")
for (i in 1:(length(roiLevels) - 1)) {
  for (j in (i + 1):length(roiLevels)) {
    roi_i <- as.character(roiLevels[i])
    roi_j <- as.character(roiLevels[j])
    identical_inputs <- identical(matrixADs[[roi_i]], matrixADs[[roi_j]])
    cat(paste0("SCC_", region, "_", roi_i, " vs SCC_", region, "_", roi_j, ": ", identical_inputs, "\n"))
  }
}

# 3. Compare SCC mean functions across ROI levels
cat("\n📐 L2 Norm differences between mean curves:\n")
for (i in 1:(length(roiLevels) - 1)) {
  for (j in (i + 1):length(roiLevels)) {
    roi_i <- as.character(roiLevels[i])
    roi_j <- as.character(roiLevels[j])
    diffNorm <- sqrt(sum((meanCurves[[roi_i]] - meanCurves[[roi_j]])^2))
    cat(paste0("Mean estimate diff between ROI ", roi_i, " and ROI ", roi_j, ": ", round(diffNorm, 5), "\n"))
  }
}

# 4. Overlay SCC mean functions visually (identical curves = perfect stability)

curve_df <- data.frame(index = seq_along(meanCurves[[1]]))
for (roi in roiLevels) {
  curve_df[[paste0("ROI_", roi)]] <- meanCurves[[as.character(roi)]]
}
curve_long <- tidyr::pivot_longer(curve_df, cols = -index,
                                  names_to = "ROI_level", values_to = "mean_estimate")

ggplot(curve_long, aes(x = index, y = mean_estimate, color = ROI_level)) +
  geom_line(size = 1.1) +
  labs(
    title = paste("SCC Mean Function Across ROI Levels —", region),
    x = "Voxel index", y = "Mean PET difference",
    color = "ROI"
  ) +
  theme_minimal(base_family = "serif") +
  scale_color_brewer(palette = "Set1")

# 5. Compare detected points spatially across ROI levels
detectedList <- list()
for (roi in roiLevels) {
  scc <- sccObjects[[as.character(roi)]]
  points <- neuroSCC::getPoints(scc)$positivePoints
  points$roi <- roi
  detectedList[[as.character(roi)]] <- points
}
allDetected <- bind_rows(detectedList)
detectionCounts <- allDetected %>% group_by(roi) %>% summarise(n_detected = n())
print(detectionCounts)

# Overlap between ROI 1 and others
cat("\n📊 Overlap of detected points with ROI 1:\n")
for (roi in roiLevels[-1]) {
  overlap <- inner_join(detectedList[["1"]], detectedList[[as.character(roi)]], by = c("x", "y"))
  cat(paste0("ROI ", roi, ": ", nrow(overlap), " shared points\n"))
}

# 6. Visualize overlap of SCC detections and ground truth ROI
roiFile <- file.path("roisNormalizadas", paste0("wwwx", region, "_redim_crop_squ_flipLR_newDim_", subjectID, ".nii"))
truePoints <- neuroSCC::processROIs(
  roiFile = roiFile,
  region = region,
  number = subjectID,
  save = FALSE,
  verbose = FALSE
)
roiPositive <- subset(truePoints, z == paramZ & pet == 1, select = c("x", "y"))

ggplot(allDetected, aes(x = x, y = y, color = factor(roi))) +
  geom_point(size = 1.6, alpha = 0.7) +
  geom_point(data = roiPositive, aes(x = x, y = y), color = "black", shape = 15, size = 2) +
  coord_fixed() +
  labs(title = paste("SCC Detected Points —", region, subjectID),
       color = "ROI Level",
       subtitle = "Black = True ROI") +
  theme_minimal(base_family = "serif")
