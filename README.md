# Simultaneous Confidence Corridors in Neuroimage Data Analysis: Applications for Alzheimer's Disease Diagnosis

This repository contains all the code and data used in our article (currently under review):

Arias-López, J. A., Cadarso-Suárez, C., & Aguiar-Fernández, P. (Under review). Simultaneous Confidence Corridors in neuroimage data analysis: applications for Alzheimer's Disease diagnosis.

The analysis combines both MATLAB and R workflows:
- MATLAB for initial neuroimaging pre-processing (realignment, unwrapping, coregistration, normalization, masking)
- R scripts for SCC computation and SPM comparison analysis

All code is documented to ensure reproducibility. Feel free to use this code to replicate our results or adapt it for your own research. Long live open science!

![Comparative analysis of detection capabilities](Auxiliary%20Files/ppv_spm_30.png)
*Comparative analysis of detection capabilities between Simultaneous Confidence Corridors (SCC, red) and Statistical Parametric Mapping (SPM, blue) across different brain regions*

The analysis of single-patient vs group comparison can be found in the continuation of this work: [PhD-2024-SCC-vs-SPM-SinglePatient-vs-Group](https://github.com/iguanamarina/PhD-2024-SCC-vs-SPM-SinglePatient-vs-Group)