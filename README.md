# PCA_projection
 Project genotypes into PCA space using pre-made SNP loadings.
 This repo is currently under construction.

## Make PC Loadings
- From reference data, extract sites in reference data that are common (MAF>0.02) present in TOPMed, and in approximate linkage equilibrium
- Remove related individuals (<3rd degree, using --king-cutoff in plink2) from refernece data
- Run PCA, save loadings

## Project Samples on PCs
- Extract sites in PC loadings, require that __% of sites from loadings be present in data or fail
- Use plink2 --score
