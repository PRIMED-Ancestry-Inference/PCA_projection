version 1.0

task subsetVariants {
	input {
		File vcf
		Array[File] variant_files = []
		File? sample_file
		File? alt_allele_file
		Float? min_maf
		Float? missingness_filter
		Int genome_build = 38
		Boolean snps_only = true
		Boolean rm_dup = true
		String output_chr = "chrM"
		Int mem_gb = 8
	}

	Int disk_size = ceil(2.5*(size(vcf, "GB") + size(variant_files, "GB"))) + 5
	String filename = basename(vcf)
	String basename = if (sub(filename, ".bcf", "") != filename) then basename(filename, ".bcf") else basename(filename, ".vcf.gz")
	String prefix = if (sub(filename, ".bcf", "") != filename) then "--bcf" else "--vcf"

	command <<<
		#get list of ranges to exclude
		wget https://raw.githubusercontent.com/GrindeLab/PCA/refs/heads/main/data/highLD/exclude_b~{genome_build}.txt
		cut -f 1,2,3 exclude_b~{genome_build}.txt > exclude.txt

		#subset file with --extract extract.txt
		plink2 \
			~{prefix} ~{vcf} \
			~{"--maf " + min_maf} \
			~{if length(variant_files) > 0 then "--extract-intersect " else ""} ~{sep=" " variant_files} \
			~{"--keep " + sample_file} \
			~{"--geno " + missingness_filter} \
			--exclude bed1 exclude.txt \
			~{true="--snps-only 'just-acgt' --max-alleles 2" false="" snps_only} \
			~{true="--rm-dup force-first" false="" rm_dup} \
			~{"--alt1-allele 'force' " + alt_allele_file + " 2 1 '#'"} \
			--allow-extra-chr \
			--chr 1-22 \
			--output-chr ~{output_chr} \
			--set-all-var-ids @:#:\$r:\$a \
			--double-id \
			--make-bed --out ~{basename}_subset
		awk '/^[^#]/ {print $2}' ~{basename}_subset.bim > selected_variants.txt
	>>>

	output {
		File snps_to_keep="selected_variants.txt"
		File subset_bed="~{basename}_subset.bed"
		File subset_bim="~{basename}_subset.bim"
		File subset_fam="~{basename}_subset.fam"
		File subset_log="~{basename}_subset.log"
	}

	runtime {
		docker: "quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
	}
}


#prune dataset by linkage
task pruneVars {
	input{
		File bed
		File bim
		File fam
		Int window_size = 10000
		Int shift_size = 1000
		Float r2_threshold = 0.1
		String output_chr = "chrM"
		Int mem_gb = 8
	}

	Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB")))
	String basename = basename(bed, ".bed")

	command <<<
		set -e -o pipefail

		command="plink2 --bed ~{bed} --bim ~{bim} --fam ~{fam} \
			--rm-dup force-first \
			--output-chr ~{output_chr} \
			--set-all-var-ids @:#:\$r:\$a \
			--indep-pairwise ~{window_size} ~{shift_size} ~{r2_threshold} \
			--out ~{basename}_indep"
		printf "${command}\n"
		${command}

		# extract pruned variants
		command="plink2 --bed ~{bed} --bim ~{bim} --fam ~{fam} \
			--extract ~{basename}_indep.prune.in \
			--output-chr ~{output_chr} \
			--set-all-var-ids @:#:\$r:\$a \
			--make-bed \
			--out ~{basename}_pruned"
		printf "${command}\n"
		${command}
	>>>

	output {
		#File subset_keep_vars="~{basename}_indep.prune.in"
		File out_bed="~{basename}_pruned.bed"
		File out_bim="~{basename}_pruned.bim"
		File out_fam="~{basename}_pruned.fam"
	}

	runtime {
		docker: "quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
	}
}
