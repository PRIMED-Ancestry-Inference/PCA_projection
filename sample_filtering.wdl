version 1.0

#remove related individuals
task removeRelateds {
	input {
		File bed
		File bim
		File fam
		Float max_kinship_coefficient = 0.0442
		Int mem_gb = 16
	}

	Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB"))) + 10
	String basename = basename(bed, ".bed")

	command <<<
		#identify individuals who are less related than kinship threshold
		command="plink2 --bed ~{bed} --bim ~{bim} --fam ~{fam} \
		--king-cutoff ~{max_kinship_coefficient} \
		--output-chr chrM \
		--set-all-var-ids @:#:\$r:\$a \
		--make-bed \
		--out ~{basename}_unrel"
		printf "${command}\n"
		${command}
	>>>

	output {
		#File subset_keep_inds="~{basename}.king.cutoff.in.id"
		File out_bed="~{basename}_unrel.bed"
		File out_bim="~{basename}_unrel.bim"
		File out_fam="~{basename}_unrel.fam"
	}

	runtime {
		docker: "quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
	}
}


# convert to numeric chromosomes before running king
task king_ibdseg {
	input {
		File bed
		File bim
		File fam
		Int degree = 3
		Int mem_gb = 16
		Int n_cpus = 4
	}

	Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB"))) + 10
	String basename = basename(bed, ".bed")

	command <<<
		set -e -o pipefail

		plink --bed ~{bed} --bim ~{bim} --fam ~{fam} \
			--output-chr 26 \
			--make-bed --out tmp \

		king -b tmp.bed \
			--ibdseg --degree ~{degree} \
			--prefix ~{basename} \
			--cpus ~{n_cpus}

		Rscript -e "\
		library(dplyr); \
		library(readr); \
		kin <- read_tsv('~{basename}.seg'); \
		kin <- mutate(kin, IBD0=(1 - IBD1Seg - IBD2Seg), Kinship=0.5*PropIBD); \
		write_tsv(kin, '~{basename}.kin0'); \
		"
	>>>

	output {
		File kin0 = "~{basename}.kin0"
	}

	runtime {
		docker: "uwgac/topmed-master:2.12.1"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
		cpu: n_cpus
	}
}

# convert to numeric chromosomes before running king
task king_robust {
	input {
		File bed
		File bim
		File fam
		Int degree = 3
		Int mem_gb = 16
		Int n_cpus = 4
	}

	Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB"))) + 10
	String basename = basename(bed, ".bed")

	command <<<
		set -e -o pipefail

		plink --bed ~{bed} --bim ~{bim} --fam ~{fam} \
			--output-chr 26 \
			--make-bed --out tmp \

		king -b tmp.bed \
			--kinship --degree ~{degree} \
			--prefix ~{basename} \
			--cpus ~{n_cpus}

	>>>

	output {
		File kin = "~{basename}.kin"
	}

	runtime {
		docker: "uwgac/topmed-master:2.12.1"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
		cpu: n_cpus
	}
}


task findRelated {
	input {
		File king_file
		String estimator = "PropIBD"
		Int degree = 3
		Int mem_gb = 16
	}

	command <<<
		Rscript -e "\
		library(GENESIS); \
		exponent <- c(5,7,9,11,13)[~{degree}]; \
		thresh <- 2^(-exponent/2); \
		kinobj <- kingToMatrix('~{king_file}', estimator='~{estimator}', thresh=thresh); \
		part <- pcairPartition(kinobj, kin.thresh=thresh); \
		rels <- part[['rels']]; \
		readr::write_tsv(tibble::tibble(FID=rels, IID=rels), 'related_samples.txt', col_names=FALSE); \
		"
	>>>

	output {
		File related_samples = "related_samples.txt"
	}

	runtime {
		docker: "uwgac/topmed-master:2.12.1"
		memory: mem_gb + " GB"
	}
}


task removeSamples {
	input {
		File bed
		File bim
		File fam
		File samples_to_remove
		String suffix = "subset"
		Int mem_gb = 16
	}

	Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB"))) + 10
	String basename = basename(bed, ".bed")

	command <<<
		command="plink2 --bed ~{bed} --bim ~{bim} --fam ~{fam} \
		--remove ~{samples_to_remove} \
		--output-chr chrM \
		--set-all-var-ids @:#:\$r:\$a \
		--make-bed \
		--out ~{basename}_~{suffix}"
		printf "${command}\n"
		${command}
	>>>

	output {
		File out_bed="~{basename}_~{suffix}.bed"
		File out_bim="~{basename}_~{suffix}.bim"
		File out_fam="~{basename}_~{suffix}.fam"
	}

	runtime {
		docker: "quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
	}
}
