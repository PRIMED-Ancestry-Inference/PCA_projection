version 1.0

task identifyColumns {
	input {
		File ref_variants
		String id_column = "ID"
	}

	command <<<
		Rscript -e "\
		dat <- readr::read_tsv('~{ref_variants}', comment = '##', n_max=100); \
		if (ncol(dat) == 1) id_col <- 1 else id_col <- which(names(dat) == '~{id_column}'); \
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
		Array[File] bed
		Array[File] bim
		Array[File] fam
		Boolean rm_dup = true
		Int mem_gb = 16
	}

	Int disk_size = ceil(3*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB"))) + 10

	command <<<
		# merge plink files
		cat ~{write_lines(bed)} | sed 's/.bed//' > bfile.txt
		plink2 --pmerge-list bfile.txt bfile \
			--merge-max-allele-ct 2 \
			--out tmp
		plink2 --bfile tmp \
			--output-chr chrM \
			--set-all-var-ids @:#:\$r:\$a \
			--make-bed --out merged \
			~{true="--rm-dup exclude-all" false="" rm_dup}
	>>>

	output {
		File out_bed = "merged.bed"
		File out_bim = "merged.bim"
		File out_fam = "merged.fam"
	}

	runtime {
		docker: "quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
	}
}
