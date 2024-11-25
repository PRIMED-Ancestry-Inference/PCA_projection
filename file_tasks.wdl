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
		docker: "quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
	}
}


task pgen2bed {
    input {
        File pgen
        File pvar
        File psam
        File? alt_allele_file
        String? out_prefix
        Int mem_gb = 16
    }

    Int disk_size = ceil(3*(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB"))) + 10
    String out_string = if defined(out_prefix) then out_prefix else basename(pgen, ".pgen")

    command {
        plink2 \
            --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} \
            ~{"--alt1-allele 'force' " + alt_allele_file + " 2 1 '#'"} \
            --make-bed \
            --out ${out_string}
        md5sum ${out_string}.bed | cut -d " " -f 1 > md5_bed.txt
        md5sum ${out_string}.bim | cut -d " " -f 1 > md5_bim.txt
        md5sum ${out_string}.fam | cut -d " " -f 1 > md5_fam.txt
    }

    output {
        File out_bed = "${out_string}.bed"
        File out_bim = "${out_string}.bim"
        File out_fam = "${out_string}.fam"
        Map[String, String] md5sum = {
            "bed": read_string("md5_bed.txt"), 
            "bim": read_string("md5_bim.txt"), 
            "fam": read_string("md5_fam.txt")
        }
    }

    runtime {
        docker: "quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_gb + " GB"
    }
}
