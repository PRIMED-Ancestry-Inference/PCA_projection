version 1.0

#remove related individuals
task removeRelateds {
	input {
		File bed
		File bim
		File fam
		Float max_kinship_coefficient = 0.0442
		Int mem_gb = 8
	}

	#Float disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB"))) * 1.5	#hoping this works?
	String basename = basename(bed, ".bed")

	command <<<
		#make the kinship matrix- #this is recommended when accessing the matrix numerous times but may be less efficient here
		#command="/plink2 --bed ~{bed} --bim ~{bim} --fam ~{fam} \
		#--make-king triangle bin \
		#--out ref_kin"
		#printf "${command}\n"
		#${command}

		#identify individuals who are less related than kinship threshold
		command="/plink2 --bed ~{bed} --bim ~{bim} --fam ~{fam} \
		--king-cutoff ~{max_kinship_coefficient} \
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
		docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
		#disks: "local-disk " + disk_size + " HDD"
		memory: mem_gb + " GB"
	}
}

#extract variants in common with some database (don't want to create PC's with SNPs that won't be found in most datasets)
task extractOverlap {
	input{
		File ref_bim
		File bed
		File bim
		File fam
		Int mem_gb = 8
	}

	#Float disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB"))) * 1.5	#hoping this works?
	String basename = basename(bed, ".bed")

	command <<<
		#get a list of variant names in common between the two, save to extract.txt
		awk 'FNR==NR{a[$2]; next}{if($2 in a){print $2}}' ~{ref_bim} ~{bim} > extract.txt

		#subset bed with --extract extract.txt
		command="/plink2 --bed ~{bed} --bim ~{bim} --fam ~{fam} \
			--keep-allele-order \
			--extract extract.txt \
			--make-bed \
			--out ~{basename}_overlap"
		printf "${command}\n"
		${command}
	>>>

	output {
		File snps_to_keep="extract.txt"
		File subset_bed="~{basename}_overlap.bed"
		File subset_bim="~{basename}_overlap.bim"
		File subset_fam="~{basename}_overlap.fam"
		File subset_log="~{basename}_overlap.log"
	}

	runtime {
		docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
		#disks: "local-disk " + disk_size + " HDD"
		memory: mem_gb + " GB"
	}
}

#prune dataset by linkage
task pruneVars {
	input{
		File bed
		File bim
		File fam
		#File keep_inds
		Int window_size = 10000
		Int shift_size = 1000
		Float r2_threshold = 0.1
		Int mem_gb = 8
	}

	#Float disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB"))) * 1.5	#hoping this works?
	String basename = basename(bed, ".bed")
	
	command <<<
		command="/plink2 --bed ~{bed} --bim ~{bim} --fam ~{fam} \
			--keep-allele-order \
			--indep-pairwise ~{window_size} ~{shift_size} ~{r2_threshold} \
			--out ~{basename}_indep"
		printf "${command}\n"
		${command}

		# extract pruned variants
		command="/plink2 --bed ~{bed} --bim ~{bim} --fam ~{fam} \
			--extract ~{basename}_indep.prune.in \
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
		docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
		#disks: "local-disk " + disk_size + " HDD"
		memory: mem_gb + " GB"
	}
}

task make_pca_loadings {
	input {
		File bed
		File bim
		File fam
		#File keep_inds
		#File keep_vars
	}

	#Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB")))
	String basename = basename(bed, ".bed")
	#ln --symbolic ${P} ${basename}.${k}.P.in

	command <<<
		command="/plink2 --bed ~{bed} --bim ~{bim} --fam ~{fam} \
			--freq counts \
			--pca allele-wts \
			--out ~{basename}_snp_loadings"
		printf "${command}\n"
		${command}
	>>>

	runtime {
		docker: "emosyne/plink2@sha256:195614c953e81da763661be20ef149be7d16b348cb68c5d54114e261aede1c92"
		#docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
		#disks: "local-disk " + disk_size + " HDD"
		#memory: mem_gb + " GB"
		#cpu: n_cpus
	}

	output {
		#check output file name from --score in plink2
		File var_freq_counts = "~{basename}_snp_loadings.acount"
		File snp_loadings = "~{basename}_snp_loadings.eigenvec.allele" 
		File projection_log = "~{basename}_snp_loadings.log"
	}
}

#run the reference panel with loadings just created- can take most of this from other workflow
task run_pca_projected {
	input {
		File bed
		File bim
		File fam
		#File keep_vars
		File loadings
		File freq_file
	}

	#Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB")))
	String basename = basename(bed, ".bed")

	command <<<
		#https://www.cog-genomics.org/plink/2.0/score#pca_project
#		command="/plink2 --bed ${basename} \
		command="/plink2 --bed ~{bed} --bim ~{bim} --fam ~{fam} \
			--read-freq ~{freq_file} \
			--score ~{loadings} 2 5 header-read no-mean-imputation variance-standardize \
			--score-col-nums 6-15 \
			--out ~{basename}_proj_pca"
		printf "${command}\n"
		${command}
	>>>

	runtime {
		docker: "emosyne/plink2@sha256:195614c953e81da763661be20ef149be7d16b348cb68c5d54114e261aede1c92"
		#docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
		#disks: "local-disk " + disk_size + " HDD"
		#memory: mem_gb + " GB"
		#cpu: n_cpus
	}

	output {
		#check output file name from --score in plink2
		File pca_projection = "~{basename}_proj_pca.sscore"
		File projection_log = "~{basename}_proj_pca.log"
	}
}

workflow create_pca_projection {
	input{ 
		File bed
		File bim
		File fam
		File ref_bim
   	 	Boolean remove_relateds = true
		Float? max_kinship_coefficient
    	Boolean prune_variants = true
		Int? window_size
		Int? shift_size
		Int? r2_threshold
		#String? mem_gb
		#Int? n_cpus
	}

  	if (remove_relateds) {
		call removeRelateds {
			input:
				bed = bed,
				bim = bim,
				fam = fam,
				max_kinship_coefficient = max_kinship_coefficient
		}
	}

	call extractOverlap {
		input:
			ref_bim = ref_bim,
			bed = select_first([removeRelateds.out_bed, bed]),
			bim = select_first([removeRelateds.out_bim, fam]),
			fam = select_first([removeRelateds.out_bim, fam]),
	}

	if (prune_variants) {
		call pruneVars {
			input:
				bed = extractOverlap.subset_bed,
				bim = extractOverlap.subset_bim,
				fam = extractOverlap.subset_fam,
				window_size = window_size,
				shift_size = shift_size,
				r2_threshold = r2_threshold
		}
	}

	File final_bed = select_first([pruneVars.out_bed, extractOverlap.subset_bed])
	File final_bim = select_first([pruneVars.out_bim, extractOverlap.subset_bim])
	File final_fam = select_first([pruneVars.out_fam, extractOverlap.subset_fam])

	call make_pca_loadings {
		input:
			bed = final_bed,
			bim = final_bim,
			fam = final_fam
	}

	call run_pca_projected {
		input:
#			bed = extractOverlap.subset_bed, #might be able to just send the original .bed file here
#			bim = extractOverlap.subset_bim,
#			fam = extractOverlap.subset_fam,
			bed = final_bed,
			bim = final_bim, 
			fam = final_fam,
			loadings = make_pca_loadings.snp_loadings,
			freq_file = make_pca_loadings.var_freq_counts,
			#mem_gb = mem_gb,
			#n_cpus = n_cpus
	}

	meta {
		author: "Jonathan Shortt"
		email: "jonathan.shortt@cuanschutz.edu"
		description: "## run_projected_pca\n This workflow is used to create a pca projection from a genetic reference dataset (in plink format, i.e., .bed/.bim/.fam). First, the reference data is subsetted to include only sites in common with a provided reference .bim (intended to contain only variants that one would expect to find in all downstream datsets that will be projected using loadings created in this worflow (e.g., a list of common sites that are easily imputed in TOPMed)), and then pruned for linkage equilibrium (after removing related individuals). Then pca is run on the dataset."
	}
}
