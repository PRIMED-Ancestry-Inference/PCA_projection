version 1.0

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
		docker: "quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0"
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
		docker: "quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
		cpu: n_cpus
	}

	output {
		File projection_file = "~{basename}_proj_pca.sscore"
		File projection_log = "~{basename}_proj_pca.log"
	}
}
