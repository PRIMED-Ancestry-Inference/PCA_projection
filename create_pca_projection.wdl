version 1.0

import "variant_filtering.wdl" as variant_tasks
import "sample_filtering.wdl" as sample_tasks
import "file_tasks.wdl" as file_tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/palantir-workflows/main/ImputationPipeline/PCATasks.wdl" as pca_tasks
import "pca_plots.wdl" as pca_plots

workflow create_pca_projection {
	input{ 
		Array[File] vcf
		File? ref_variants
		File? sample_file
		Float? missingness_filter
		Int? n_pcs
		Int? genome_build
		Boolean prune_variants = true
		Boolean remove_relateds = true
		Float? min_maf
		Float? max_kinship_coefficient
		Int? window_size
		Int? shift_size
		Float? r2_threshold
		File? groups_file
	}

	if (defined(ref_variants)) {
		call file_tasks.identifyColumns {
			input:
				ref_variants = select_first([ref_variants, ""])
		}
	}

	scatter (file in vcf) {
		call variant_tasks.subsetVariants {
			input:
				vcf = file,
				variant_file = identifyColumns.id_file,
				sample_file = sample_file,
				missingness_filter = missingness_filter,
				genome_build = genome_build,
				min_maf = min_maf
		}

		if (prune_variants) {
			call variant_tasks.pruneVars {
				input:
					pgen = subsetVariants.subset_pgen,
					pvar = subsetVariants.subset_pvar,
					psam = subsetVariants.subset_psam,
					window_size = window_size,
					shift_size = shift_size,
					r2_threshold = r2_threshold
			}
		}

		File subset_pgen = select_first([pruneVars.out_pgen, subsetVariants.subset_pgen])
		File subset_pvar = select_first([pruneVars.out_pvar, subsetVariants.subset_pvar])
		File subset_psam = select_first([pruneVars.out_psam, subsetVariants.subset_psam])
	}

	if (length(vcf) > 1) {
		call file_tasks.mergeFiles {
			input:
				pgen = subset_pgen,
				pvar = subset_pvar,
				psam = subset_psam
		}
	}

	File merged_pgen = select_first([mergeFiles.out_pgen, pruneVars.out_pgen[0], subsetVariants.subset_pgen[0]])
	File merged_pvar = select_first([mergeFiles.out_pvar, pruneVars.out_pvar[0], subsetVariants.subset_pvar[0]])
	File merged_psam = select_first([mergeFiles.out_psam, pruneVars.out_psam[0], subsetVariants.subset_psam[0]])

  	if (remove_relateds) {
		call sample_tasks.removeRelateds {
			input:
				pgen = merged_pgen,
				pvar = merged_pvar,
				psam = merged_psam,
				max_kinship_coefficient = max_kinship_coefficient
		}
	}

	File final_pgen = select_first([removeRelateds.out_pgen, merged_pgen])
	File final_pvar = select_first([removeRelateds.out_pvar, merged_pvar])
	File final_psam = select_first([removeRelateds.out_psam, merged_psam])	

	call file_tasks.pgen2bed {
		input:
			pgen = final_pgen,
			pvar = final_pvar,
			psam = final_psam
	}

	call pca_tasks.PerformPCA {
		input:
			bed = pgen2bed.out_bed,
			bim = pgen2bed.out_bim,
			fam = pgen2bed.out_fam,
			basename = basename(pgen2bed.out_bed, ".bed"),
			n_pcs = n_pcs
	}

	call pca_plots.run_pca_plots {
		input: 
			data_file = PerformPCA.pcs, 
			groups_file = groups_file
	}

	output {
		File pcs = PerformPCA.pcs
		File pc_variance = PerformPCA.pc_variance
		File pc_loadings = PerformPCA.pc_loadings
		File mean_sd = PerformPCA.mean_sd
		File eigenvectors = PerformPCA.eigenvectors
		File eigenvalues = PerformPCA.eigenvalues
		File? pca_plots_pc12 = run_pca_plots.pca_plots_pc12
		Array[File]? pca_plots_pairs = run_pca_plots.pca_plots_pairs
		File? pca_plots_parcoord = run_pca_plots.pca_plots_parcoord
		File? pca_plots = run_pca_plots.pca_plots
	}

	meta {
		author: "Jonathan Shortt, Stephanie Gogarten, Amy Watt"
		email: "jonathan.shortt@cuanschutz.edu"
		description: "This workflow is used to create a pca projection from a genetic reference dataset (in VCF format). First, the reference data is subsetted to include only sites in common with a provided reference variant file (intended to contain only variants that one would expect to find in all downstream datsets that will be projected using loadings created in this worflow (e.g., a list of common sites that are easily imputed in TOPMed)), and then pruned for linkage equilibrium. The related individuals are removed. Then PCA is run on the dataset."
	}
}
