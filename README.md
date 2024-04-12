# PCA_projection

Project genotypes into PCA space using pre-made SNP loadings.

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


Outputs:

output | description
--- | ---
var_freq_counts | counts of variant frequencies
snp_loadings | SNP loadings
loadings_log | log from running plink2 --pca
pca_projection | PCs from running PCA on this dataset with calculated loadings
projection_log | log from running plink2 --score


## projected_PCA

This workflow is used to project a genetic test dataset (in VCF format) into PCA space using user-defined allele loadings. First, the allele loadings (from the create_pca_projection workflow) and the test dataset are both subset to contain the same set of variants (Note: this workflow assumes that variants from both the loadings and test dataset have been previously harmonized such that variants follow the same naming convention.) Then the test dataset is projected onto the principal components.


Inputs:

input | description
--- | ---
ref_loadings | File with SNP loadings (e.g. snp_loadings output from create_pca_projection)
ref_freqs | File with variant frequencies (e.g. var_freq_counts output from create_pca_projection)
vcf | Array of VCF files (possibly split by chromosome)
min_overlap | minimum overlap between variants in loadings and vcf files (default 0.95). If the overlap is less than this threshold, PCA will not be run and the workflow will exit.


Outputs:

output | description
--- | ---
projection_file | PCs from running PCA on this dataset with ref_loadings
projection_log | log from running plink2 --score



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
