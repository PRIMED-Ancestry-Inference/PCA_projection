version 1.0

import "variant_filtering.wdl" as variant_tasks
import "file_tasks.wdl" as file_tasks
import "pca_tasks.wdl" as pca_tasks

workflow projected_PCA {
	input {
		File ref_loadings
		File ref_freqs
		Array[File] vcf
		Float min_overlap = 0.95
	}

	call identifyColumns {
		input:
			ref_loadings = ref_loadings
	}

	scatter (file in vcf) {
		call variant_tasks.subsetVariants {
			input:
				vcf = file,
				variant_file = ref_loadings,
				variant_id_col = identifyColumns.id_col
		}
	}

	if (length(vcf) > 1) {
		call file_tasks.mergeFiles {
			input:
				pgen = subsetVariants.subset_pgen,
				pvar = subsetVariants.subset_pvar,
				psam = subsetVariants.subset_psam
		}
	}

	File final_pgen = select_first([mergeFiles.out_pgen, subsetVariants.subset_pgen[0]])
	File final_pvar = select_first([mergeFiles.out_pvar, subsetVariants.subset_pvar[0]])
	File final_psam = select_first([mergeFiles.out_psam, subsetVariants.subset_psam[0]])

	call checkOverlap {
		input:
			variant_file = ref_loadings,
			pvar = final_pvar
	}

	#check for overlap, if overlap is less than threshold, stop
	if (checkOverlap.overlap >= min_overlap) {
		call pca_tasks.run_pca_projected {
			input:
				pgen = final_pgen,
				pvar = final_pvar,
				psam = final_psam,
				loadings = ref_loadings,
				freq_file = ref_freqs,
				id_col = identifyColumns.id_col,
				allele_col = identifyColumns.allele_col,
				pc_col_first = identifyColumns.pc_col_first,
				pc_col_last = identifyColumns.pc_col_last
		}
	}

	output {
		File? projection_file = run_pca_projected.projection_file
		File? projection_log = run_pca_projected.projection_log
		Float overlap = checkOverlap.overlap
	}

	meta {
		author: "Jonathan Shortt, Stephanie Gogarten"
		email: "jonathan.shortt@cuanschutz.edu"
		description: "This workflow is used to project a genetic test dataset (in VCF format) into PCA space using user-defined allele loadings. First, the allele loadings (from the create_pca_projection workflow) and the test dataset are both subset to contain the same set of variants (Note: this workflow assumes that variants from both the loadings and test dataset have been previously harmonized such that variants follow the same naming convention, alleles at each site are ordered identically, and variants are sorted). Then the test dataset is projected onto the principal components."
	}
}


task identifyColumns {
	input {
		File ref_loadings
	}

	command <<<
		Rscript -e "\
		dat <- readr::read_tsv('~{ref_loadings}', n_max=100)
		writeLines(as.character(which(names(dat) == 'ID')), 'id_col.txt')
		if (is.element('A1', names(dat))) allele_col <- 'A1' else allele_col <- 'ALT'
		writeLines(as.character(which(names(dat) == allele_col)), 'allele_col.txt')
		pc_cols <- grep('^PC', names(dat))
		writeLines(as.character(pc_cols[1]), 'pc_col_first.txt')
		writeLines(as.character(pc_cols[length(pc_cols)]), 'pc_col_last.txt')
		"
	>>>

	output {
		Int id_col = read_int("id_col.txt")
		Int allele_col = read_int("allele_col.txt")
		Int pc_col_first = read_int("pc_col_first.txt")
		Int pc_col_last = read_int("pc_col_last.txt")
	}

	runtime {
		docker: "us.gcr.io/broad-dsp-gcr-public/anvil-rstudio-bioconductor:3.17.0"
	}
}


task checkOverlap {
	input {
		File variant_file
		File pvar
	}

	command <<<
		python3 <<CODE
		def countLines (myfile):
			line_count=0
			with open(myfile, "r") as infp:
				for line in infp:
					if not (line.startswith("#")):
						line_count=line_count + 1
			return line_count
		
		loadings_count=countLines("~{variant_file}")
		new_loadings_count=countLines("~{pvar}")
		proportion=float(new_loadings_count)/loadings_count
		print("%.3f" % proportion)
		CODE
	>>>

	output {
		Float overlap = read_float(stdout())
	}

	runtime {
		docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
	}
}
