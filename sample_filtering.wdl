version 1.0

#remove related individuals
task removeRelateds {
	input {
		File pgen
		File pvar
		File psam
		File? king_table
		Float max_kinship_coefficient = 0.0442
		Int mem_gb = 16
	}

	Int disk_size = ceil(1.5*(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB"))) + 10
	String basename = basename(pgen, ".pgen")

	command <<<
		#identify individuals who are less related than kinship threshold
		command="plink2 --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} \
		--king-cutoff ~{max_kinship_coefficient} ~{'---king-cutoff-table ' + king_table} \
		--output-chr chrM \
		--set-all-var-ids @:#:\$r:\$a \
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
		docker: "quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
	}
}


task king {
	input {
		File bed
		File bim
		File fam
		Int degree = 4
		Int mem_gb = 16
		Int n_cpus = 4
	}

	Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB"))) + 10
	String basename = basename(bed, ".bed")

	command <<<
		king -b ~{bed} \
			--related --degree ~{degree} \
			--prefix ~{basename}_unrel \
			--cpus ~{n_cpus}
	>>>

	output {
		File kin0 = "~{basename}_unrel.kin0"
	}

	runtime {
		docker: "uwgac/topmed-master:2.12.1"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
		cpu: n_cpus
	}
}
