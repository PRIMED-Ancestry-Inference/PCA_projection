version 1.0

task subsetVariants {
	input {
		File vcf
		File? variant_file
		Float? min_maf
		Int genome_build = 38
		Boolean snps_only = true
		Boolean rm_dup = true
		Boolean set_var_ids = true
		Int mem_gb = 8
	}

	Int disk_size = ceil(2.5*(size(vcf, "GB"))) + 5
	String filename = basename(vcf)
	String basename = if (sub(filename, ".bcf", "") != filename) then basename(filename, ".bcf") else basename(filename, ".vcf.gz")
	String prefix = if (sub(filename, ".bcf", "") != filename) then "--bcf" else "--vcf"

	command <<<
		#get list of ranges to exclude
		wget https://raw.githubusercontent.com/GrindeLab/PCA/refs/heads/main/data/highLD/exclude_b~{genome_build}.txt
		cut -f 1,2,3 exclude_b~{genome_build}.txt > exclude.txt

		#subset file with --extract extract.txt
		plink2 ~{prefix} ~{vcf} ~{"--maf " + min_maf} ~{"--extract " + variant_file} \
			--exclude bed1 exclude.txt \
			~{true="--snps-only 'just-acgt'" false="" snps_only} \
			~{true="--rm-dup force-first" false="" rm_dup} \
			--output-chr chrM \
			~{true="--set-all-var-ids @:#:\$r:\$a" false="" set_var_ids} \
			--make-pgen --out ~{basename}_subset
		awk '/^[^#]/ {print $3}' ~{basename}_subset.pvar > selected_variants.txt
	>>>

	output {
		File snps_to_keep="selected_variants.txt"
		File subset_pgen="~{basename}_subset.pgen"
		File subset_pvar="~{basename}_subset.pvar"
		File subset_psam="~{basename}_subset.psam"
		File subset_log="~{basename}_subset.log"
	}

	runtime {
		docker: "quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
	}
}


#prune dataset by linkage
task pruneVars {
	input{
		File pgen
		File pvar
		File psam
		Int window_size = 10000
		Int shift_size = 1000
		Float r2_threshold = 0.1
		Int mem_gb = 8
	}

	Int disk_size = ceil(1.5*(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB")))
	String basename = basename(pgen, ".pgen")
	
	command <<<
		command="plink2 --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} \
			--rm-dup force-first \
			--output-chr chrM \
			--set-all-var-ids @:#:\$r:\$a \
			--indep-pairwise ~{window_size} ~{shift_size} ~{r2_threshold} \
			--out ~{basename}_indep"
		printf "${command}\n"
		${command}

		# extract pruned variants
		command="plink2 --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} \
			--extract ~{basename}_indep.prune.in \
			--output-chr chrM \
			--set-all-var-ids @:#:\$r:\$a \
			--make-pgen \
			--out ~{basename}_pruned"
		printf "${command}\n"
		${command}
	>>>

	output {
		#File subset_keep_vars="~{basename}_indep.prune.in"
		File out_pgen="~{basename}_pruned.pgen"
		File out_pvar="~{basename}_pruned.pvar"
		File out_psam="~{basename}_pruned.psam"
	}

	runtime {
		docker: "quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
	}
}
