version 1.0

workflow projected_PCA {
	input {
		File ref_loadings
		File ref_freqs
		File vcf
		Float overlap = 0.95
		Int? mem_gb
		Int? n_cpus
	}

	call prepareFiles {
		input:
			ref_loadings = ref_loadings,
			ref_freqs = ref_freqs,
			vcf = vcf
		}

	call checkOverlap {
		input:
			ref_loadings = ref_loadings,
			pca_loadings = prepareFiles.subset_loadings
	}

	if (checkOverlap.overlap >= overlap) {
		call run_pca_projected {
			input:
				pgen = prepareFiles.subset_pgen,
				pvar = prepareFiles.subset_pvar,
				pfam = prepareFiles.subset_pfam,
				loadings = prepareFiles.subset_loadings,
				freq_file = prepareFiles.subset_freqs,
				mem_gb = mem_gb,
				n_cpus = n_cpus
		}
	}

	output {
		File? projection_file = run_pca_projected.projection_file
		File? projection_log = run_pca_projected.projection_log
	}

	meta {
		author: "Jonathan Shortt"
		email: "jonathan.shortt@cuanschutz.edu"
		description: "This workflow is used to project a genetic test dataset (in plink format, i.e., .bed/.bim/.fam) into PCA space using user-defined allele loadings. First, the allele loadings (from the create_pca_projection workflow) and the test dataset are both subset to contain the same set of variants (Note: this workflow assumes that variants from both the loadings and test dataset have been previously harmonized such that variants follow the same naming convention, alleles at each site are ordered identically, and variants are sorted). Then the test dataset is projected onto the principal components."
	}
}

task prepareFiles {
	input {
		File ref_loadings
		File ref_freqs
		File vcf
		Int mem_gb = 8
	}

	Int disk_size = ceil(2.5*(size(vcf, "GB")))
	String filename = basename(vcf)
	String b1 = sub(filename, "bcf", "")
	String b2 = sub(filename, ".bcf", "")
	String v1 = sub(filename, "vcf.gz", "")
	String v2 = sub(filename, ".vcf.gz", "")
	String basename = if (sub(filename, ".bcf", "") != filename) then basename(filename, ".bcf") else basename(filename, ".vcf.gz")
	String in_file = if (sub(filename, ".bcf", "") != filename) then "--bcf " + vcf else "--vcf " + vcf

	command <<<
		#get a list of variant names in common between the two, save to extract.txt
		#variant name in loadings is assumed to be 3rd column, assuming plink2 format (https://www.cog-genomics.org/plink/2.0/formats#eigenvec)
		awk '{print $2}' ~{ref_loadings} > extract.txt
		#subset file with --extract extract.txt
		/plink2 ~{in_file} --extract extract.txt --make-pgen --out ~{basename}_pcaReady

		#extract variants in-common variants from ref_loadings
		#this step may not be necessary at all since plink --score might just be able to deal with it
		head -n 1 ~{ref_loadings} > loadings_pcaReady.txt
		awk 'FNR==NR{a[$1]; next}{if($2 in a) {print $0}}' extract.txt ~{ref_loadings} >> loadings_pcaReady.txt

		#extract variants in-common variants from ref_freqs
		#this step may not be necessary at all since plink --score might just be able to deal with it
		head -n 1 ~{ref_freqs} > freqs_pcaReady.txt
		awk 'FNR==NR{a[$1]; next}{if($2 in a) {print $0}}' extract.txt ~{ref_freqs} >> freqs_pcaReady.txt
	>>>

	output {
		File snps_to_keep="extract.txt"
		File subset_pgen="~{basename}_pcaReady.pgen"
		File subset_pvar="~{basename}_pcaReady.pvar"
		File subset_pfam="~{basename}_pcaReady.pfam"
		File subset_log="~{basename}_pcaReady.log"
		File subset_loadings="loadings_pcaReady.txt"
		File subset_freqs="freqs_pcaReady.txt"
	}

	runtime {
		docker: "emosyne/plink2@sha256:195614c953e81da763661be20ef149be7d16b348cb68c5d54114e261aede1c92"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
	}
}

task checkOverlap {
	input {
		File ref_loadings
		File pca_loadings
		Int mem_gb = 8
	}

	command <<<
		python3 <<CODE
		def countLines (myfile):
			line_count=0
			with open(myfile, "r") as infp:
				for line in infp:
					line_count=line_count + 1
			return line_count
		
		loadings_count=countLines("~{ref_loadings}")
		new_loadings_count=countLines("~{pca_loadings}")
		proportion=float(loadings_count)/new_loadings_count
		print("%.3f" % proportion)
		CODE
	>>>

	#check for overlap, if overlap is less than threshold, stop, default overlap threshold is 0.95

	output {
		Float overlap = read_float(stdout())
	}

	runtime {
		docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
		memory: mem_gb + " GB"
	}
}

task run_pca_projected {
	input {
		File pgen
		File pvar
		File pfam
		File loadings
		File freq_file
		Int mem_gb = 8
		Int n_cpus = 4 # check this
	}

	Int disk_size = ceil(1.5*(size(pgen, "GB") + size(pvar, "GB") + size(pfam, "GB")))
	String basename = basename(pgen, ".pgen")
	#ln --symbolic ${P} ${basename}.${k}.P.in

	command <<<
		#https://www.cog-genomics.org/plink/2.0/score#pca_project
		command="/plink2 --pgen ~{pgen} --pvar ~{pvar} --pfam ~{pfam} \
			--read-freq ~{freq_file} \
			--score ~{loadings} 2 5 header-read no-mean-imputation variance-standardize \
			--score-col-nums 6-15 \
			--out ~{basename}_proj_pca"
		printf "${command}\n"
		${command}
	>>>

	runtime {
		docker: "emosyne/plink2@sha256:195614c953e81da763661be20ef149be7d16b348cb68c5d54114e261aede1c92"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
		cpu: n_cpus
	}

	output {
		#check output file name from --score in plink2
		File projection_file = "~{basename}_proj_pca.sscore"
		File projection_log = "~{basename}_proj_pca.log"
	}
}
