version 1.0

task prepareFiles {
    input {
        File ref_loadings
        File ref_freqs
        File bed
        File bim
        File fam
        Float? overlap
        Int mem_gb = 8
    }

    Float disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB"))) * 1.5    #hoping this works?
    String basename = basename(bed, ".bed")
    
    command <<<
        #get a list of variant names in common between the two, save to extract.txt
        #variant name in loadings is assumed to be 3rd column, assuming plink2 format (https://www.cog-genomics.org/plink/2.0/formats#eigenvec)
        awk 'FNR==NR{a[$3]; next}{if($2 in a){print $2}}' ~{ref_loadings} ~{bim} > extract.txt

        #subset bed with --extract extract.txt
        /plink2 --bed ~{bed} --bim ~{bim} --fam ~{fam} --extract extract.txt --make-bed --keep-allele-order --out ~{basename}_pcaReady
        
        #extract variants in-common variants from ref_loadings
        #this step may not be necessary at all since plink --score might just be able to deal with it
        head -n 1 ~{ref_loadings} > loadings_pcaReady.txt
        awk 'FNR==NR{a[$1]; next}{if($3 in a) {print $0}}' extract.txt ~{ref_loadings} >> loadings_pcaReady.txt
    
        #extract variants in-common variants from ref_freqs
        #this step may not be necessary at all since plink --score might just be able to deal with it
        head -n 1 ~{ref_freqs} > freqs_pcaReady.txt
        awk 'FNR==NR{a[$1]; next}{if($3 in a) {print $0}}' extract.txt ~{ref_freqs} >> freqs_pcaReady.txt
    
        #check for overlap, if overlap is less than threshold, stop, default overlap threshold is 0.95
        loadings_count=$(tail -n +2  < ${ref_loadings} | wc -l)
        new_loadings_count=$(tail -n +2 < loadings_pcaReady.txt | wc -l)
        #doing this calculation in awk because I'm lousy at bash...
        prop=$(awk -v old=${loadings_count} -v new=${new_loadings_count} '{prop=new/old; print prop}' )
        printf "Variant overlap is ${prop} of original.\n"
        #https://support.terra.bio/hc/en-us/articles/360037484851-Variable-Types-in-WDL#:~:text=When%20working%20with%20optional%20variables%20in%20your%20command%2C,The%20syntax%20for%20that%20is%3A%20%24%20%7Bdefault%3D%22value%22%20variableName%7D
        myoverlap=${default=0.95 overlap}
        exit_code=$(awk -v prop=${prop} -v default_threshold=${myoverlap} '{myexit=0; if(prop < default_threshold){myexit=1}; print myexit }' )
        if [${exit_code} -gt 0]; then
            printf "SNP overlap ${prop} is lower than" # fix sentence fragment here
            exit ${exit_code};
        fi
    >>>
    
    output {
        File snps_to_keep="extract.txt"
        File subset_bed="~{basename}_pcaReady.bed"
        File subset_bim="~{basename}_pcaReady.bim"
        File subset_fam="~{basename}_pcaReady.fam"
        File subset_log="~{basename}_pcaReady.log"
        File subset_loadings="loadings_pcaReady.txt"
        File subset_freqs="freqs_pcaReady.txt"
    }
    
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
        #disks: "local-disk " + disk_size + " HDD"
        memory: mem_gb + " GB"
    }
}

task run_pca_projected {
    input {
        File bed
        File bim
        File fam
        File loadings
        File freq_file
        Int mem_gb = 8
        Int n_cpus = 4 # check this
      }

    Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB")))
    String basename = basename(bed, ".bed")
    #ln --symbolic ${P} ${basename}.${k}.P.in

    command <<<
        #https://www.cog-genomics.org/plink/2.0/score#pca_project
        command="/plink2 --bfile ${basename} \
            --read-freq ~{freq_file} \
            --score ~{loadings} 2 5 header-read no-mean-imputation variance-standardize \
                   --score-col-nums 6-15 \
            --out ${basename}_proj_pca"
        printf "${command}\n"
        ${command}
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
        #disks: "local-disk " + disk_size + " HDD"
        memory: mem_gb + " GB"
        cpu: n_cpus
      }

    output {
        #check output file name from --score in plink2
        File pca_projection = "~{basename}_pca.sscore"
        File projection_log = "~{basename}_pca.log"
    }
}

workflow pca_projection {
    input {
        File ref_loadings
        File ref_freqs
        File bed
        File bim
        File fam
        Int? mem_gb
        Int? n_cpus
        
    }
    
    call prepareFiles {
        input:
            ref_loadings = ref_loadings,
            ref_freqs = ref_freqs,
            bed = bed, 
            bim = bim, 
            fam = fam
    }
    
    call run_pca_projected {
        input:
            bed = prepareFiles.subset_bed,
                bim = prepareFiles.subset_bim,
                fam = prepareFiles.subset_fam,
                loadings = prepareFiles.subset_loadings,
            freq_file = prepareFiles.subset_freqs,
                mem_gb = mem_gb,
                n_cpus = n_cpus
    }
    
    meta {
            author: "Jonathan Shortt"
            email: "jonathan.shortt@cuanschutz.edu"
            description: "## run_projected_pca\n This workflow is used to project a genetic test dataset (in plink format, i.e., .bed/.bim/.fam) into pca space using user-defined allele loadings. First, the allele loadings (.P produced by ADMIXTURE) and the test dataset are both subset to contain the same set of variants (Note: this workflow assumes that variants from both the loadings and test dataset have been previously harmonized such that variants follow the same naming convention, alleles at each site are ordered identically, and variants are sorted). Then the test dataset is projected onto the principal components."
    }

}
