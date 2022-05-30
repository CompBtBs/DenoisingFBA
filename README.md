# Combining denoising of RNA-seq data and Flux Balance Analysis for cluster analysis of single cells

 All RNA-seq data are available from the Gene Expression Omnibus database (accession numbers GSE110949, GSE118056) or  EBI Single Cell Expression Atlas (accession number E-GEOD-86618).
 
 The datasets must be downloaded as saved in the directory data in H5 format (from Scanpy).
 
 The output directory contains:
 -the Reaction Activity Scores values for all the datasets and the denoisers
 -the metabolic fluxes computering using Flux Balance Analysis (FBA)
 
 To run the codes you need to install:
 * Scanpy
 * Cobrapy
 * MAGIC denoiser
 * ENHANCE scripts
 * SAVER (R-based)
 
 We provided an R script to run the SAVER denoiser.
 
 The notebook "workflow can be used to compute both RAS from transcriptional data and single cell Flux Balance Analysis
