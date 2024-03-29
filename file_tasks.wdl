version 1.0

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
