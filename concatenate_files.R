library(argparser)
library(readr)
library(tidyr)

# Get parameters 
argp <- arg_parser("colormap")
argp <- add_argument(parser = argp, 
                     arg = "--ref_pcs",
                     type = "character", 
                     nargs = 1, 
                     help="Tab-delimited file with IID column and PCs for reference data")
argp <- add_argument(parser = argp, 
                     arg = "--ref_groups",
                     type = "character", 
                     nargs = 1, 
                     help="Two-column tab-delimited file with subject_id column and group label")
argp <- add_argument(parser = argp, 
                     arg = "--projection_file",
                     type = "character", 
                     nargs = 1, 
                     help="Tab-delimited file with IID column and PCs for sample data")

argv <- parse_args(argp)
ref_pcs <- argv$ref_pcs
ref_groups <- argv$ref_groups
projection_file <- argv$projection_file

# Read data
projection_file <- as.data.frame(read_tsv(projection_file))

if(is.na(ref_pcs)) {
  merged_pcs <- projection_file
  merged_groups <- data.frame(subject_id = projection_file$IID, group = rep("projected_samples", nrow(projection_file)))
} else {
  
  # Read ref data 
  ref_pcs <- as.data.frame(read_tsv(ref_pcs))
  
  # Cast IID as character
  ref_pcs$IID <- as.character(ref_pcs$IID)
  projection_file$IID <- as.character(projection_file$IID)
  
  # Concatenate files
  merged_pcs <- rbind(ref_pcs, projection_file)
  
  # Define colors 
  paired <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
  dark2 <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666')
  
  if(!is.na(ref_groups)) {
    ref_groups <- as.data.frame(read_tsv(ref_groups))
    ref_groups$subject_id <- as.character(ref_groups$subject_id)
    groups <- unique(ref_groups$group)
    n_groups <- length(groups)
    
    if(n_groups <= 7) {
      palette <- dark2[1:n_groups]
    } else if(n_groups <= 12) {
      palette <- paired[1:ln_groups]
    } else {
      ggplotColours <- function(n = 6, h = c(0, 360) + 15){
        if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
        hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
      }
      
      palette <- ggplotColours(n=n_groups)
    } 
    
    colormap <- tibble(
      group = c(groups, "projected_samples"),
      color = c(palette[1:n_groups], '#666666')
    )
    df <- data.frame(subject_id = projection_file$IID, group = rep("projected_samples", nrow(projection_file)))
    colnames(df) <- c("subject_id", "group")
    
    merged_groups <- rbind(ref_groups, df)
    
  } else {
    colormap <- tibble(
      group = c("ref_samples", "projected_samples"),
      color = c('#1f78b4', '#666666')
    )
    
    merged_groups <- rbind(data.frame(subject_id = ref_pcs$IID, group = rep("ref_samples", nrow(ref_pcs))), 
                           data.frame(subject_id = projection_file$IID, group = rep("projected_samples", nrow(projection_file))), 
                           by = "subject_id")
    colnames(merged_groups) <- c("subject_id", "group")
  }
}


write_tsv(colormap, "colormap.tsv")
write_tsv(merged_pcs, "merged_pcs.tsv")
write_tsv(merged_groups, "merged_groups.tsv")




# Rscript concatenate_files.R --ref_pcs test_data/ref_pcs.tsv --ref_groups test_data/ref_groups.tsv --projection_file test_data/pca_plots_test_data.sscore
# Rscript concatenate_files.R --ref_pcs ref_proj_pca.sscore --ref_groups 1000G_populations.txt --projection_file proj_pca.sscore