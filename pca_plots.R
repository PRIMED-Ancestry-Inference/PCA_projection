library(readr)
library(ggplot2)
library(GGally)
library(tidyr)
library(argparser)

# Rscript pca_plots.R --data_file test_data/pca_plots_test_data.sscore --n_pairs 3

# Get parameters 
argp <- arg_parser("PCA plots")
argp <- add_argument(parser = argp, 
                     arg = "--data_file",
                     type = "character", 
                     nargs = 1, 
                     help="path to data file")

argp <- add_argument(parser = argp, 
                     arg = "--group",
                     type = "character", 
                     nargs = 1, 
                     help="Variable name for grouping results")

argp <- add_argument(parser = argp, 
                     arg = "--n_pairs",
                     type = "integer", 
                     nargs = 1, 
                     help="Number of PCs to use for pairwise plots")

argv <- parse_args(argp)
data_file <- argv$data_file
group <- argv$group
n_pairs <- argv$n_pairs

# Read data file 
dat <- as.data.frame(read_tsv(data_file))
colnames(dat)[colnames(dat)=="#IID"] <- "IID"
rownames(dat) <- dat[,1]
dat <- dat[,-(1:3)]
pc_columns <- grep("^PC", names(dat), value = TRUE)
num_pcs <- length(pc_columns)

print(num_pcs)

# Color by group 
if(!is.na(argv$group)) {
  # join to table by sample id to find groups 
  # stopifnot group %in% colnames(table_w_groups)
  # dat <- left_join(dat, table_w_groups, by = "IID")
} else {
  group <- "group"
  dat$group <- "NA"
}

# pca <- pivot_longer(dat, cols = c(starts_with("PC")), names_to = "pc", values_to = "value")

# PC 1 and PC 2
p <- ggplot(dat, aes(PC1_AVG, PC2_AVG, color=group)) +
  geom_point(alpha=0.5) +
  theme(legend.position="none")
ggsave("test_data/out_file_pc12.pdf", plot = p, width = 7, height = 6)

# Pairwise PCs from inputs
npr <- min(num_pcs, n_pairs)
p <- ggpairs(dat, mapping = aes(color = .data[[group]], alpha = 0.5), columns = 1:npr)
ggsave("test_data/out_file_pairs.pdf", plot = p, width=8, height=8)

# Parcoord Plots
p <- ggparcoord(dat, columns = 1:num_pcs, groupColumn = group, alphaLines = 0.5, scale = "uniminmax") +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=2))) +
  xlab("PC") + ylab("")
ggsave("test_data/out_file_parcoord.pdf", plot = p, width = 7, height = 6)
