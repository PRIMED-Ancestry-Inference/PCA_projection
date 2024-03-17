version 1.0

import "projected_pca.wdl" as tasks

workflow create_pca_projection {
	input{ 
		Array[File] vcf
		File ref_variants
    	Boolean prune_variants = true
   	 	Boolean remove_relateds = true
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
		call tasks.subsetVariants {
			input:
				vcf = file,
				variant_file = ref_variants,
				variant_id_col = identifyColumns.id_col
		}

		if (prune_variants) {
			call pruneVars {
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
		call tasks.mergeFiles {
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
		call removeRelateds {
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

	call make_pca_loadings {
		input:
			pgen = final_pgen,
			pvar = final_pvar,
			psam = final_psam
	}

	call tasks.run_pca_projected {
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
		dat <- readr::read_tsv('~{ref_variants}', comment = '##', n_max=100)
		writeLines(as.character(which(names(dat) == 'ID')), 'id_col.txt')
		"
	>>>

	output {
		Int id_col = read_int("id_col.txt")
	}

	runtime {
		docker: "us.gcr.io/broad-dsp-gcr-public/anvil-rstudio-bioconductor:3.17.0"
	}
}


#remove related individuals
task removeRelateds {
	input {
		File pgen
		File pvar
		File psam
		Float max_kinship_coefficient = 0.0442
		Int mem_gb = 8
	}

	Int disk_size = ceil(1.5*(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB")))
	String basename = basename(pgen, ".pgen")

	command <<<
		#identify individuals who are less related than kinship threshold
		command="/plink2 --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} \
		--king-cutoff ~{max_kinship_coefficient} \
		--make-pgen \
		--out ~{basename}_unrel"
		printf "${command}\n"
		${command}
	>>>

	output {
		#File subset_keep_inds="~{basename}.king.cutoff.in.id"
		File out_pgen="~{basename}_unrel.pgen"
		File out_pvar="~{basename}_unrel.pvar"
		File out_psam="~{basename}_unrel.psam"
	}

	runtime {
		docker: "emosyne/plink2@sha256:195614c953e81da763661be20ef149be7d16b348cb68c5d54114e261aede1c92"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
	}
}


#prune dataset by linkage
task pruneVars {
	input{
		File pgen
		File pvar
		File psam
		Int window_size = 10000
		Int shift_size = 1000
		Float r2_threshold = 0.1
		Int mem_gb = 8
	}

	Int disk_size = ceil(1.5*(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB")))
	String basename = basename(pgen, ".pgen")
	
	command <<<
		command="/plink2 --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} \
			--rm-dup force-first --set-missing-var-ids @:#:\$r:\$a \
			--indep-pairwise ~{window_size} ~{shift_size} ~{r2_threshold} \
			--out ~{basename}_indep"
		printf "${command}\n"
		${command}

		# extract pruned variants
		command="/plink2 --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} \
			--extract ~{basename}_indep.prune.in \
			--make-pgen \
			--out ~{basename}_pruned"
		printf "${command}\n"
		${command}
	>>>

	output {
		#File subset_keep_vars="~{basename}_indep.prune.in"
		File out_pgen="~{basename}_pruned.pgen"
		File out_pvar="~{basename}_pruned.pvar"
		File out_psam="~{basename}_pruned.psam"
	}

	runtime {
		docker: "emosyne/plink2@sha256:195614c953e81da763661be20ef149be7d16b348cb68c5d54114e261aede1c92"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
	}
}


task make_pca_loadings {
	input {
		File pgen
		File pvar
		File psam
		Int mem_gb = 8
		Int n_cpus = 4
	}

	Int disk_size = ceil(1.5*(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB")))
	String basename = basename(pgen, ".pgen")

	command <<<
		command="/plink2 --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} \
			--freq counts \
			--pca allele-wts \
			--out ~{basename}_snp_loadings"
		printf "${command}\n"
		${command}
	>>>

	runtime {
		docker: "emosyne/plink2@sha256:195614c953e81da763661be20ef149be7d16b348cb68c5d54114e261aede1c92"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
		cpu: n_cpus
	}

	output {
		File var_freq_counts = "~{basename}_snp_loadings.acount"
		File snp_loadings = "~{basename}_snp_loadings.eigenvec.allele" 
		File projection_log = "~{basename}_snp_loadings.log"
	}
}
