# mapggm_supplemental: Supplemental files for Multi-attribute perturbed Gaussian graphical models 

This repository contains files used in ``Detection of multiple perturbations in multi-omics biological networks'' (Griffin, Johnson, & Kolaczyk; under revision).

## Pipeline
The **pipeline** directory includes active R scripts, job management, postprocessing, and plotting files for the simulation study described in the paper.

The R package **mapggm** called from these files contains functions for multi-attribute network estimation and perturbation detection.  The source for this package can be found at <https://github.com/paulajgriffin/mapggm>.  To use this package, install the devtools package from CRAN and run: 

	library(devtools)
	install_github('paulajgriffin/mapggm')

## TCGA data
The **TCGA** directory contains preprocessing code and data files used in Section 5 of the main paper.


