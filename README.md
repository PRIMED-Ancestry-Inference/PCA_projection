# PCA_projection

Run PCA and project genotypes into PCA space using pre-made SNP loadings.

## create_pca_projection

This workflow is used to create a pca projection from a genetic reference dataset (in VCF format). First, the reference data is subsetted to include only sites in common with a provided reference variant file (intended to contain only variants that one would expect to find in all downstream datsets that will be projected using loadings created in this worflow (e.g., a list of common sites that are easily imputed in TOPMed)), and then pruned for linkage equilibrium. The related individuals are removed. Then PCA is run on the dataset.

Inputs:

input | description
--- | ---
vcf | Array of VCF files (possibly split by chromosome)
ref_variants | file with variants to use in the PCA calculation. The column with variant IDs should be labeled 'ID'.
prune_variants | Boolean for whether to do LD pruning on the variants (default true)
min_maf | minimum MAF for variants to include (optional)
remove_relateds | Boolean for whether to remove samples with relatedness above max_kinship_coefficient (default true)
max_kinship_coefficient | if remove_relateds is true, remove one of each pair of samples with kinship > this value (default 0.0442 = 3rd degree relatives)
window_size | window size for LD pruning (default 10,000)
shift_size | shift size for LD pruning (default 1000)
r2_threshold | r2 threshold for LD pruning (default 0.1)
groups_file | Two-column tsv file of subject_id and group, used to label plots (optional)


Outputs:

output | description
--- | ---
pcs | PCs for samples used to create projection
pc_variance | variance explained by each PC
pc_loadings | SNP loadings
mean_sd | mean and SD for each variant in SNP loadings file
eigenvectors | eigenvectors
eigenvalues | eigenvalues
loadings_log | log from running plink2 --pca
pca_projection | PCs from running PCA on this dataset with calculated loadings
projection_log | log from running plink2 --score
pca_plots_pc12 | png file of PC1 and PC2 scatterplot
pca_plot_pairs | png file of pairwise PC scatterplots 
pca_plots_parcoord | png file of parallel coordinates plot for PCs
pca_plots | html file with PCA plots 


## projected_PCA

This workflow is used to project a genetic test dataset (in VCF format) into PCA space using user-defined allele loadings. First, the allele loadings (from the create_pca_projection workflow) and the test dataset are both subset to contain the same set of variants (Note: this workflow assumes that variants from both the loadings and test dataset have been previously harmonized such that variants follow the same naming convention.) Then the test dataset is projected onto the principal components.


Inputs:

input | description
--- | ---
ref_loadings | File with SNP loadings (e.g. pc_loadings output from create_pca_projection)
ref_meansd | File with variant mean and SD (e.g. mean_sd output from create_pca_projection)
ref_pcs | PCs from running PCA on reference dataset to create joint plots (optional)
ref_groups | Two-column tsv file of subject_id and group from reference dataset, used to label plots (optional)
groups_file | Two-column tsv file of subject_id and group from sample dataset, used to label plots (optional)
vcf | Array of VCF files (possibly split by chromosome)
min_overlap | minimum overlap between variants in loadings and vcf files (default 0.95). If the overlap is less than this threshold, PCA will not be run and the workflow will exit.


Outputs:

output | description
--- | ---
projection_file | PCs from running PCA on this dataset with ref_loadings
projection_log | log from running plink2 --score
pca_plots_pc12 | png file of PC1 and PC2 scatterplot of samples
pca_plot_pairs | png file of pairwise PC scatterplots of samples
pca_plots_parcoord | png file of parallel coordinates plot for PCs of samples
pca_plots | html file with PCA plots of samples 
pca_plots_pc12_ref | png file of PC1 and PC2 scatterplot of samples overlaid on references 
pca_plot_pairs_ref | png file of pairwise PC scatterplots of samples overlaid on references 
pca_plots_parcoord_ref | png file of parallel coordinates plot for PCs of samples overlaid on references 
pca_plots_ref | html file with PCA plots of samples overlaid on references 


## LD_pruning

This workflow prunes variants for linkage equilibrium.

Inputs:

input | description
--- | ---
vcf | Array of VCF files (possibly split by chromosome)
variant_file | Optional file with variant selection to start the pruning
variant_id_col | Column in variant_file containing the IDs
min_maf | minimum MAF for variants to include (optional)
snps_only | Boolean for whether to use only SNPs (default true) 
window_size | window size for LD pruning (default 10,000)
shift_size | shift size for LD pruning (default 1000)
r2_threshold | r2 threshold for LD pruning (default 0.1)

Outputs:

output | description
--- | ---
pruned_vcf | Array of pruned VCF files



## select_variants_by_pop_maf

This workflow selects all variants with MAF > a minimum threshold in any population (i.e. the union of filtering by MAF in each population separately). Samples to select for each population are identified by reading the population_descriptor and sample tables from the specified workspace.

The output of this workflow is a text file with variant IDs, taken from the ID column of the VCF file(s). Any missing values in the ID column are replaced with chr:pos:ref:alt. Duplicate variant IDs are excluded.


Inputs:

input | description
--- | ---
vcf | Array of VCF files (possibly split by chromosome)
min_maf | minimum MAF for variants to select
population_descriptor | the descriptor to use for identifying populations
population_labels | Array of labels for each population. If this input is not supplied, the workflow will use all unique labels for the population descriptor.
workspace_name | name of the workspace with a population_descriptor data table (e.g. "PRIMED_1000G")
workspace_namespace | namespace of the workspace (e.g. "primed-data-cc-1")


Outputs:

output | description
--- | ---
maf_filtered_variants | Text file with variants that passed the MAF filter in any population.



## pca_plots

This workflow is uses a file with PCs to create pairs plots and parallel coordinate plots. 

Inputs:

input | description
--- | ---
data_file | PCs from running PCA (ie pcs from create_pca_projection or projection_file from projected_pca)
groups_file | Two-column tsv file of subject_id and group, used to label plots (optional)
colormap | Two-column tsv file of group and color, used to color plots (optional)
n_pairs | number of PCs to use for pairs plots (default 10)

Outputs:

output | description
--- | ---
pca_plots_pc12 | png file of PC1 and PC2 scatterplot
pca_plot_pairs | png file of pairwise PC scatterplots 
pca_plots_parcoord | png file of parallel coordinates plot for PCs
pca_plots | html file with PCA plots 
