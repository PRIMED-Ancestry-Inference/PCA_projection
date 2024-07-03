library(argparser)

# Rscript pca_plots.R --data_file test_data/pca_plots_test_data.sscore --groups_file test_data/groups_file_test.tsv --n_pairs 3 --path_to_rmd ~/Downloads/PCA_projection 
# Rscript pca_plots.R --data_file test_data/pca_plots_test_data.sscore --n_pairs 3 --path_to_rmd ~/Downloads/PCA_projection 
# Rscript pca_plots.R --data_file merged_pcs.tsv --groups_file merged_groups.tsv --colormap colormap.tsv --n_pairs 3 --path_to_rmd ~/Downloads/PCA_projection 

# Rscript pca_plots.R --data_file test_data/pca_plots_test_data.sscore --groups_file test_data/groups_file_test.tsv --n_pairs 3 --path_to_rmd /home/rstudio/PCA_projection
# Rscript pca_plots.R --data_file test_data/pca_plots_test_data.sscore $(if [ -f "test_data/groups_file_test.tsv" ]; then echo "--groups_file test_data/groups_file_test.tsv"; fi) --n_pairs 3 --path_to_rmd ~/Downloads/PCA_projection


# Get parameters 
argp <- arg_parser("PCA plots")
argp <- add_argument(parser = argp, 
                     arg = "--data_file",
                     type = "character", 
                     nargs = 1, 
                     help="Path to data file")

argp <- add_argument(parser = argp, 
                     arg = "--groups_file",
                     type = "character", 
                     nargs = 1, 
                     help="Two-column tab-delimited file with subject ID and group label")

argp <- add_argument(parser = argp, 
                     arg = "--colormap",
                     type = "character", 
                     nargs = 1, 
                     help="Two-column tab-delimited file with group and color (one row per group)")

argp <- add_argument(parser = argp, 
                     arg = "--n_pairs",
                     type = "integer", 
                     nargs = 1, 
                     help="Number of PCs to use for pairwise plots")

argp <- add_argument(parser = argp, 
                     arg = "--path_to_rmd", 
                     type = "character", 
                     nargs = 1,
                     help="rmd filepath")

argv <- parse_args(argp)
data_file <- argv$data_file
groups_file <- argv$groups_file
colormap <- argv$colormap
n_pairs <- argv$n_pairs
path_to_rmd <- argv$path_to_rmd

parameters <- list(data_file=data_file, groups_file=groups_file, colormap=colormap, n_pairs=n_pairs)

file.copy(file.path(path_to_rmd, "pca_plots.Rmd"), "pca_plots.Rmd")
rmarkdown::render(input = "pca_plots.Rmd", params = parameters, quiet=TRUE)