version 1.0

import "variant_filtering.wdl" as variant_tasks
import "sample_filtering.wdl" as sample_tasks
import "file_tasks.wdl" as file_tasks
import "pca_tasks.wdl" as pca_tasks

workflow create_pca_projection {
	input{ 
		Array[File] vcf
		File ref_variants
		Boolean prune_variants = true
		Boolean remove_relateds = true
		Float? min_maf
		Float? max_kinship_coefficient
		Int? window_size
		Int? shift_size
		Int? r2_threshold
	}

	call identifyColumns {
		input:
			ref_variants = ref_variants
	}

	scatter (file in vcf) {
		call variant_tasks.subsetVariants {
			input:
				vcf = file,
				variant_file = ref_variants,
				variant_id_col = identifyColumns.id_col,
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

	call pca_tasks.make_pca_loadings {
		input:
			pgen = final_pgen,
			pvar = final_pvar,
			psam = final_psam
	}

	call pca_tasks.run_pca_projected {
		input:
			pgen = merged_pgen,
			pvar = merged_pvar,
			psam = merged_psam,
			loadings = make_pca_loadings.snp_loadings,
			freq_file = make_pca_loadings.var_freq_counts,
			id_col = 2,
			allele_col = 5,
			pc_col_first = 6,
			pc_col_last = 15
	}

	output {
		File var_freq_counts = make_pca_loadings.var_freq_counts
		File snp_loadings =  make_pca_loadings.snp_loadings
		File loadings_log =  make_pca_loadings.projection_log
		File pca_projection = run_pca_projected.projection_file
		File projection_log = run_pca_projected.projection_log
	}

	meta {
		author: "Jonathan Shortt, Stephanie Gogarten"
		email: "jonathan.shortt@cuanschutz.edu"
		description: "This workflow is used to create a pca projection from a genetic reference dataset (in VCF format). First, the reference data is subsetted to include only sites in common with a provided reference variant file (intended to contain only variants that one would expect to find in all downstream datsets that will be projected using loadings created in this worflow (e.g., a list of common sites that are easily imputed in TOPMed)), and then pruned for linkage equilibrium. The related individuals are removed. Then PCA is run on the dataset."
	}
}


task identifyColumns {
	input {
		File ref_variants
	}

	command <<<
		Rscript -e "\
		dat <- readr::read_tsv('~{ref_variants}', comment = '##', n_max=100); \
		if (ncol(dat) == 1) id_col <- 1 else id_col <- which(names(dat) == 'ID'); \
		writeLines(as.character(id_col), 'id_col.txt')
		"
	>>>

	output {
		Int id_col = read_int("id_col.txt")
	}

	runtime {
		docker: "us.gcr.io/broad-dsp-gcr-public/anvil-rstudio-bioconductor:3.17.0"
	}
}
