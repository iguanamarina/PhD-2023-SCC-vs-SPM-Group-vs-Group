#' generateCompareRTriplets
#'
#' Generate voxel-wise binary classifications for SCC, SPM, and ROI ground truth,
#' aligned to a shared (x, y) voxel grid at a specific axial slice (e.g., z = 35).
#' Intended as input to the `compareR()` function from the testCompareR package.
#'
#' @param grid A data.frame with columns x and y, listing all voxel positions in the slice.
#'             Typically generated using `expand.grid()` from the mask dimensions.
#' @param sccCoords A data.frame of SCC-detected positive voxel coordinates with columns x and y.
#' @param spmCoords A data.frame of SPM-detected positive voxel coordinates with columns x and y.
#' @param roiCoords A data.frame of ground-truth ROI-positive voxel coordinates with columns x and y.
#'
#' @return A data.frame with 9919 rows (or as many as in grid), and three binary columns:
#'         - SCC_pred: 1 if SCC detects this voxel, 0 otherwise
#'         - SPM_pred: 1 if SPM detects this voxel, 0 otherwise
#'         - ROI_truth: 1 if voxel is truly positive in the ROI, 0 otherwise
#'
#' @examples
#' triplets <- generateCompareRTriplets(grid, sccCoords, spmCoords, roiCoords)
generateCompareRTriplets <- function(grid, sccCoords, spmCoords, roiCoords) {
  # Input checks
  stopifnot(all(c("x", "y") %in% colnames(grid)))
  stopifnot(all(c("x", "y") %in% colnames(sccCoords)))
  stopifnot(all(c("x", "y") %in% colnames(spmCoords)))
  stopifnot(all(c("x", "y") %in% colnames(roiCoords)))
  
  # Generate voxel IDs for fast matching
  voxelIDs <- interaction(grid$x, grid$y)
  sccIDs   <- interaction(sccCoords$x, sccCoords$y)
  spmIDs   <- interaction(spmCoords$x, spmCoords$y)
  roiIDs   <- interaction(roiCoords$x, roiCoords$y)
  
  # Match voxel IDs to each method's detections
  grid$SCC_pred   <- as.integer(voxelIDs %in% sccIDs)
  grid$SPM_pred   <- as.integer(voxelIDs %in% spmIDs)
  grid$ROI_truth  <- as.integer(voxelIDs %in% roiIDs)
  
  # Return only the required binary columns
  return(grid[, c("SCC_pred", "SPM_pred", "ROI_truth")])
}
