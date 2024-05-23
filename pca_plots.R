library(readr)
library(ggplot2)
library(GGally)
library(tidyr)
library(dplyr)
library(argparser)

# Rscript pca_plots.R --data_file test_data/pca_plots_test_data.sscore --groups_file test_data/groups_file_test.tsv --n_pairs 3

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
                     arg = "--n_pairs",
                     type = "integer", 
                     nargs = 1, 
                     help="Number of PCs to use for pairwise plots")

argv <- parse_args(argp)
data_file <- argv$data_file
groups_file <- argv$groups_file
n_pairs <- argv$n_pairs

# Read data file 
dat <- as.data.frame(read_tsv(data_file))
colnames(dat)[colnames(dat)=="#IID"] <- "IID"
# rownames(dat) <- dat[,1]
# dat <- dat[,-(1:3)]
pc_columns <- grep("^PC", names(dat), value = TRUE)
num_pcs <- length(pc_columns)

# Color by group 
if(!is.na(groups_file)) {
  groups_dat <- as.data.frame(read_tsv(groups_file, col_names=TRUE, col_types="cc"))
  colnames(groups_dat)[1] <- "IID"
  group <- colnames(groups_dat)[2]
  dat <- merge(dat, groups_dat, by = "IID", all.x = TRUE)
} else {
  group <- "group"
  dat$group <- "NA"
}

# Filter data for columns relevant for PC plots 
dat <- dat %>% select(
  starts_with("PC"), 
  group
)

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
