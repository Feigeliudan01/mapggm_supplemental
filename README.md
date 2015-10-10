# mapggm_supplemental: Supplemental files for Multi-attribute perturbed Gaussian graphical models 

This repository contains files used in ``Detection of multiple perturbations in multi-omics biological networks'' (Griffin, Johnson, & Kolaczyk; under review).

## Simulations
The **simulations** directory includes active R scripts, job management, postprocessing, and plotting files for the simulation study described in the paper.

The R package **mapggm** called from these files contains functions for multi-attribute network estimation and perturbation detection.  The source for this package can be found at <https://github.com/paulajgriffin/mapggm>.  To use this package, install the devtools package from CRAN and run: 

	library(devtools)
	install_github('paulajgriffin/mapggm')

This directory contains files to be run in the following order:
- **0_DesignSims.R**: Set up the parameters for simulations (number of replicates, partial correlations, network density, etc)
- **1_SimulateCovariances.R**: Simulate networks and control data based on the defined parameters
- **2_EstimatePrecisions.R**: Estimate networks (precision/covariance) based on the control data 
- **3_EvaluateSims.R**: Generate case data and evalute performance of estimation and testing procedure
- **4_PlotSims.R**: Plot the results and generate LaTeX tables 

## TCGA data
The **TCGA** directory contains preprocessing code and data files used in Section 5 of the main paper.
The files within data_processed/Y_* are in a format that can be used with the simulation files with very small edits to change data imports.

