version 1.0

#remove related individuals
task removeRelateds {
	input {
		File pgen
		File pvar
		File psam
		Float max_kinship_coefficient = 0.0442
		Int mem_gb = 16
	}

	Int disk_size = ceil(1.5*(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB"))) + 10
	String basename = basename(pgen, ".pgen")

	command <<<
		#identify individuals who are less related than kinship threshold
		command="plink2 --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} \
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
		docker: "quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
	}
}
