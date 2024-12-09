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
		Int kinship_degree_filter = 3
		Int? window_size
		Int? shift_size
		Float? r2_threshold
		File? groups_file
		String relatedness_estimator = "robust"
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
					bed = subsetVariants.subset_bed,
					bim = subsetVariants.subset_bim,
					fam = subsetVariants.subset_fam,
					window_size = window_size,
					shift_size = shift_size,
					r2_threshold = r2_threshold
			}
		}

		File subset_bed = select_first([pruneVars.out_bed, subsetVariants.subset_bed])
		File subset_bim = select_first([pruneVars.out_bim, subsetVariants.subset_bim])
		File subset_fam = select_first([pruneVars.out_fam, subsetVariants.subset_fam])
	}

	if (length(vcf) > 1) {
		call file_tasks.mergeFiles {
			input:
				bed = subset_bed,
				bim = subset_bim,
				fam = subset_fam
		}
	}

	File merged_bed = select_first([mergeFiles.out_bed, pruneVars.out_bed[0], subsetVariants.subset_bed[0]])
	File merged_bim = select_first([mergeFiles.out_bim, pruneVars.out_bim[0], subsetVariants.subset_bim[0]])
	File merged_fam = select_first([mergeFiles.out_fam, pruneVars.out_fam[0], subsetVariants.subset_fam[0]])

  	if (remove_relateds) {

		# Map king relatedness estimator to the appropriate value for the GENESIS task.
		Map[String, String] relatedness_estimator_map = {"robust": "Kinship", "ibdseg": "PropIBD"}

		if (relatedness_estimator == "robust") {
			call sample_tasks.king_robust {
					input:
						bed = merged_bed,
						bim = merged_bim,
						fam = merged_fam,
						degree = kinship_degree_filter
				}
		}

		if (relatedness_estimator == "ibdseg") {
			call sample_tasks.king_ibdseg {
					input:
						bed = merged_bed,
						bim = merged_bim,
						fam = merged_fam,
						degree = kinship_degree_filter
				}
		}

		call sample_tasks.findRelated {
			input:
				king_file = select_first([king_robust.kin0, king_ibdseg.kin0]),
				estimator = relatedness_estimator_map[relatedness_estimator],
				degree = kinship_degree_filter
		}

		call sample_tasks.removeSamples {
			input:
				bed = merged_bed,
				bim = merged_bim,
				fam = merged_fam,
				samples_to_remove = findRelated.related_samples,
				suffix = "unrel"
		}
	}

	File final_bed = select_first([removeSamples.out_bed, merged_bed])
	File final_bim = select_first([removeSamples.out_bim, merged_bim])
	File final_fam = select_first([removeSamples.out_fam, merged_fam])

	call pca_tasks.PerformPCA {
		input:
			bed = final_bed,
			bim = final_bim,
			fam = final_fam,
			basename = basename(final_bed, ".bed"),
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
