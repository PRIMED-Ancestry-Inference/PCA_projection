version 1.0

task identifyColumns {
	input {
		File ref_variants
	}

	command <<<
		Rscript -e "\
		dat <- readr::read_tsv('~{ref_variants}', comment = '##', n_max=100); \
		if (ncol(dat) == 1) id_col <- 1 else id_col <- which(names(dat) == 'ID'); \
		system(paste('cut -f', id_col, '~{ref_variants} > variant_ids.txt')); \
		"
	>>>

	output {
		File id_file = "variant_ids.txt"
	}

	runtime {
		docker: "us.gcr.io/broad-dsp-gcr-public/anvil-rstudio-bioconductor:3.17.0"
	}
}


task mergeFiles {
	input {
		Array[File] pgen
		Array[File] pvar
		Array[File] psam
		Boolean rm_dup = true
		Int mem_gb = 16
	}

	Int disk_size = ceil(3*(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB"))) + 10

	command <<<
		# merge plink files
		cat ~{write_lines(pgen)} | sed 's/.pgen//' > pfile.txt
		plink2 --pmerge-list pfile.txt \
			--merge-max-allele-ct 2 \
			--out tmp
		plink2 --pfile tmp \
			--output-chr chrM \
			--set-all-var-ids @:#:\$r:\$a \
			--make-pgen --out merged \
			~{true="--rm-dup exclude-all" false="" rm_dup}
	>>>

	output {
		File out_pgen = "merged.pgen"
		File out_pvar = "merged.pvar"
		File out_psam = "merged.psam"
	}

	runtime {
		docker: "quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
	}
}
