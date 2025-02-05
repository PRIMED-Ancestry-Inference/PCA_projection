version 1.0

import "variant_filtering.wdl" as variant_tasks
import "file_tasks.wdl" as file_tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/palantir-workflows/main/ImputationPipeline/PCATasks.wdl" as pca_tasks
import "pca_plots.wdl" as pca_plots

workflow projected_PCA {
	input {
		File ref_loadings
		File ref_meansd
		File? ref_pcs
		File? ref_groups
		File? sample_file
		File? variant_file
		Float? missingness_filter
		File? groups_file
		Array[File] vcf
		Int? genome_build
		Float min_overlap = 0.95
	}

	call file_tasks.identifyColumns {
		input:
			ref_variants = ref_loadings,
			id_column = "SNP"
	}

	scatter (file in vcf) {
		call variant_tasks.subsetVariants {
			input:
				vcf = file,
				variant_files = select_all([identifyColumns.id_file, variant_file]),
				sample_file = sample_file,
				missingness_filter = missingness_filter,
				genome_build = genome_build
		}
	}

	if (length(vcf) > 1) {
		call file_tasks.mergeFiles {
			input:
				bed = subsetVariants.subset_bed,
				bim = subsetVariants.subset_bim,
				fam = subsetVariants.subset_fam
		}
	}

	File final_bed = select_first([mergeFiles.out_bed, subsetVariants.subset_bed[0]])
	File final_bim = select_first([mergeFiles.out_bim, subsetVariants.subset_bim[0]])
	File final_fam = select_first([mergeFiles.out_fam, subsetVariants.subset_fam[0]])

	#check for overlap, if overlap is less than threshold, stop
	call checkOverlap {
		input:
			ref_loadings = ref_loadings,
			ref_meansd = ref_meansd,
			bim = final_bim,
			min_overlap = min_overlap
	}

	call pca_tasks.ProjectArray {
		input:
			bed = final_bed,
			bim = final_bim,
			fam = final_fam,
			pc_loadings = checkOverlap.subset_loadings,
			pc_meansd = checkOverlap.subset_meansd,
			basename = basename(final_bed, ".bed")
	}

	call pca_plots.run_pca_plots {
		input:
			data_file = ProjectArray.projections,
			groups_file = groups_file
	}

	# If ref_pcs is provided, run concatenateFiles task then rerun the plotting script with output
	if (defined(ref_pcs)) {

		# need this because ref_pcs is optional but input to concatenateFiles is required
		File ref_pcs1 = select_first([ref_pcs, ""])

		call concatenateFiles {
			input:
				ref_pcs = ref_pcs1,
				ref_groups = ref_groups,
				projection_file = ProjectArray.projections
		}

		call pca_plots.run_pca_plots as run_pca_plots_ref {
			input:
				data_file = concatenateFiles.merged_pcs,
				groups_file = concatenateFiles.merged_groups,
				colormap = concatenateFiles.colormap
		}
	}

	output {
		File projection_file = ProjectArray.projections
		Float overlap = checkOverlap.overlap
		File pca_plots_pc12 = run_pca_plots.pca_plots_pc12
		Array[File] pca_plots_pairs = run_pca_plots.pca_plots_pairs
		File pca_plots_parcoord = run_pca_plots.pca_plots_parcoord
		File pca_plots = run_pca_plots.pca_plots
		File? pca_plots_pc12_ref = run_pca_plots_ref.pca_plots_pc12
		Array[File]? pca_plots_pairs_ref = run_pca_plots_ref.pca_plots_pairs
		File? pca_plots_parcoord_ref = run_pca_plots_ref.pca_plots_parcoord
		File? pca_plots_ref = run_pca_plots_ref.pca_plots
	}

	meta {
		author: "Jonathan Shortt, Stephanie Gogarten, Amy Watt"
		email: "jonathan.shortt@cuanschutz.edu"
		description: "This workflow is used to project a genetic test dataset (in VCF format) into PCA space using user-defined allele loadings. First, the allele loadings (from the create_pca_projection workflow) and the test dataset are both subset to contain the same set of variants (Note: this workflow assumes that variants from both the loadings and test dataset have been previously harmonized such that variants follow the same naming convention, alleles at each site are ordered identically, and variants are sorted). Then the test dataset is projected onto the principal components."
	}
}


task checkOverlap {
	input {
		File ref_loadings
		File ref_meansd
		File bim
		Float min_overlap
		Int mem_gb = 8
	}

	command <<<
	Rscript -e "\
	library(dplyr); \
	library(readr); \
	bim <- read_tsv('~{bim}', col_types='-c----', col_names='SNP'); \
	loadings <- read_tsv('~{ref_loadings}'); \
	new_loadings <- inner_join(bim, loadings); \
	write_tsv(new_loadings, 'subset_loadings.txt'); \
	meansd <- read_tsv('~{ref_meansd}'); \
	new_meansd <- inner_join(bim, meansd); \
	write_tsv(new_meansd, 'subset_meansd.txt'); \
	proportion <- nrow(new_loadings) / nrow(loadings); \
	prop_string <- format(proportion, digits=3); \
	writeLines(prop_string, 'overlap.txt'); \
	min_prop <- ~{min_overlap}; \
	if (proportion < min_prop) stop(paste('Variant overlap of', prop_string, 'is less than minimum of', min_prop)); \
	"
	>>>

	output {
		Float overlap = read_float("overlap.txt")
		File subset_loadings = "subset_loadings.txt"
		File subset_meansd = "subset_meansd.txt"
	}

	runtime {
		docker: "us.gcr.io/broad-dsp-gcr-public/anvil-rstudio-bioconductor:3.17.0"
		memory: mem_gb + " GB"
	}
}


task concatenateFiles {
	input {
		File ref_pcs
		File? ref_groups
		File projection_file
	}

	command <<<
	Rscript /usr/local/PCA_projection/concatenate_files.R \
		--ref_pcs ~{ref_pcs} \
		$(if [ -n "~{ref_groups}" ]; then echo "--ref_groups ~{ref_groups}"; fi) \
		--projection_file ~{projection_file}
	>>>

	output {
		File merged_pcs = "merged_pcs.tsv"
		File merged_groups = "merged_groups.tsv"
		File colormap = "colormap.tsv"
	}

	runtime{
		 docker: "uwgac/pca_projection:0.2.0"
	}
}
