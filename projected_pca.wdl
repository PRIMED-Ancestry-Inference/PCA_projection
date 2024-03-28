version 1.0

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
		call subsetVariants {
			input:
				vcf = file,
				variant_file = ref_loadings,
				variant_id_col = identifyColumns.id_col
		}
	}

	if (length(vcf) > 1) {
		call mergeFiles {
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
		call run_pca_projected {
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


task subsetVariants {
	input {
		File vcf
		File? variant_file
		Int? variant_id_col
		Boolean snps_only = true
		Boolean rm_dup = true
		Int mem_gb = 8
	}

	Int disk_size = ceil(2.5*(size(vcf, "GB"))) + 5
	String filename = basename(vcf)
	String basename = if (sub(filename, ".bcf", "") != filename) then basename(filename, ".bcf") else basename(filename, ".vcf.gz")
	String prefix = if (sub(filename, ".bcf", "") != filename) then "--bcf" else "--vcf"

	command <<<
		#get a list of variant names in common between the two, save to extract.txt
		if [ -f ~{variant_file} ]; then cut -f ~{variant_id_col} ~{variant_file} > extract.txt; fi
		#subset file with --extract extract.txt
		echo /plink2 ~{prefix} ~{vcf} ~{true="--extract extract.txt" false="" defined(variant_file)} \
			~{true="--snps-only 'just-acgt'" false="" snps_only} \
			~{true="--rm-dup force-first" false="" rm_dup} \
			--set-missing-var-ids @:#:\$r:\$a \
			--make-pgen --out ~{basename}_subset
		#awk '/^[^#]/ {print $3}' ~{basename}_subset.pvar > selected_variants.txt
	>>>

	output {
		File snps_to_keep="selected_variants.txt"
		File subset_pgen="~{basename}_subset.pgen"
		File subset_pvar="~{basename}_subset.pvar"
		File subset_psam="~{basename}_subset.psam"
		File subset_log="~{basename}_subset.log"
	}

	runtime {
		docker: "emosyne/plink2@sha256:195614c953e81da763661be20ef149be7d16b348cb68c5d54114e261aede1c92"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
	}
}


task mergeFiles {
	input {
		Array[File] pgen
		Array[File] pvar
		Array[File] psam
		Int mem_gb = 16
	}

	Int disk_size = ceil(3*(size(pgen, "GB"))) + 10

	command <<<
		# merge plink files
		cat ~{write_lines(pgen)} | sed 's/.pgen//' > pfile.txt
		/plink2 --pmerge-list pfile.txt --merge-max-allele-ct 2 --out merged
	>>>

	output {
		File out_pgen = "merged.pgen"
		File out_pvar = "merged.pvar"
		File out_psam = "merged.psam"
	}

	runtime {
		docker: "emosyne/plink2@sha256:195614c953e81da763661be20ef149be7d16b348cb68c5d54114e261aede1c92"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
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

task run_pca_projected {
	input {
		File pgen
		File pvar
		File psam
		File loadings
		File freq_file
		Int id_col
		Int allele_col
		Int pc_col_first
		Int pc_col_last
		Int mem_gb = 8
		Int n_cpus = 4
	}

	Int disk_size = ceil(1.5*(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB")))
	String basename = basename(pgen, ".pgen")

	command <<<
		#https://www.cog-genomics.org/plink/2.0/score#pca_project
		command="/plink2 --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} \
			--read-freq ~{freq_file} \
			--score ~{loadings} ~{id_col} ~{allele_col} header-read no-mean-imputation variance-standardize \
			--score-col-nums ~{pc_col_first}-~{pc_col_last} \
			--out ~{basename}_proj_pca"
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
		File projection_file = "~{basename}_proj_pca.sscore"
		File projection_log = "~{basename}_proj_pca.log"
	}
}
