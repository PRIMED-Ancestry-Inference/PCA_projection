library(argparser)
library(readr)
library(tidyr)

# Rscript colormap.R --ref_groups test_data/groups_file_test.tsv

# Get parameters 
argp <- arg_parser("colormap")
argp <- add_argument(parser = argp, 
                     arg = "--ref_groups",
                     type = "character", 
                     nargs = 1, 
                     help="Two-column tab-delimited file with subject ID and group label")

argv <- parse_args(argp)
ref_groups <- argv$ref_groups

paired <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
dark2 <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666')

if(!is.na(ref_groups)) {
  ref_groups <- as.data.frame(read_tsv(ref_groups))
  groups <- unique(ref_groups$group)
  n_groups <- length(groups)
  
  if(n_groups <= 7) {
    palette <- dark2[1:n_groups]
  } else if(n_groups <= 12) {
    palette <- paired[1:ln_groups]
  } 
  
  colormap <- tibble(
    group = c(groups, "projected_samples"),
    color = c(palette[1:n_groups], '#666666')
    )
} else {
  colormap <- tibble(
    group = c("ref_samples", "projected_samples"),
    color = c('#1f78b4', '#666666')
    )
}

write_tsv(colormap, "colormap.tsv")

