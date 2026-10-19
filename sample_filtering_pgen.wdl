version 1.0

# PLINK2 pgen/pvar/psam version of sample_filtering.wdl.
# Relatedness has to be estimated genome-wide, so per-chromosome arrays are
# first combined with mergePgen. Typical chain:
#   mergePgen -> pgen2bed -> king_ibdseg/king_robust -> findRelated -> removeSamples
#   (or mergePgen -> removeRelateds)


# merge per-chromosome filesets into one genome-wide fileset
task mergePgen {
    input {
        Array[File] pgen
        Array[File] pvar
        Array[File] psam
        String out_prefix = "merged"
        Int mem_gb = 16
    }

    Int disk_size = ceil(3*(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB"))) + 10

    command <<<
        set -e -o pipefail

        # one line per chromosome: pgen pvar psam
        paste -d ' ' ~{write_lines(pgen)} ~{write_lines(pvar)} ~{write_lines(psam)} > merge_list.txt
        cat merge_list.txt

        plink2 \
            --memory ~{mem_gb * 512} \
            --pmerge-list merge_list.txt \
            --output-chr chrM \
            --set-all-var-ids @:#:\$r:\$a \
            --make-pgen \
            --out ~{out_prefix}
    >>>

    output {
        File out_pgen = "~{out_prefix}.pgen"
        File out_pvar = "~{out_prefix}.pvar"
        File out_psam = "~{out_prefix}.psam"
    }

    runtime {
        docker: "quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_gb + " GB"
    }
}


# convert to bed with numeric chromosomes for KING
task pgen2bed {
    input {
        File pgen
        File pvar
        File psam
        String output_chr = "26"
        Int mem_gb = 16
    }

    Int disk_size = ceil(2.5*(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB"))) + 10
    String basename = basename(pgen, ".pgen")

    command <<<
        set -e -o pipefail

        plink2 --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} \
            --memory ~{mem_gb * 512} \
            --output-chr ~{output_chr} \
            --make-bed \
            --out ~{basename}
    >>>

    output {
        File out_bed = "~{basename}.bed"
        File out_bim = "~{basename}.bim"
        File out_fam = "~{basename}.fam"
    }

    runtime {
        docker: "quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_gb + " GB"
    }
}


#remove related individuals
task removeRelateds {
    input {
        File pgen
        File pvar
        File psam
        Float max_kinship_coefficient = 0.0442
        String output_chr = "chrM"
        Int mem_gb = 16
    }

    Int disk_size = ceil(1.5*(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB"))) + 10
    String basename = basename(pgen, ".pgen")

    command <<<
        set -e -o pipefail

        #identify individuals who are less related than kinship threshold
        plink2 --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} \
            --memory ~{mem_gb * 512} \
            --king-cutoff ~{max_kinship_coefficient} \
            --output-chr ~{output_chr} \
            --set-all-var-ids @:#:\$r:\$a \
            --make-pgen \
            --out ~{basename}_unrel
    >>>

    output {
        File out_pgen = "~{basename}_unrel.pgen"
        File out_pvar = "~{basename}_unrel.pvar"
        File out_psam = "~{basename}_unrel.psam"
    }

    runtime {
        docker: "quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_gb + " GB"
    }
}


# KING tasks are unchanged: give them the output of pgen2bed
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
            --memory ~{mem_gb *  512} \
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
        if (nrow(kin) > 0) kin <- mutate(kin, IBD0=(1 - IBD1Seg - IBD2Seg), Kinship=0.5*PropIBD); \
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
        File kin0 = "~{basename}.kin0"
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
        String estimator = "Kinship"
        Int degree = 3
        Int mem_gb = 16
    }

    command <<<
        Rscript -e "\
        library(GENESIS); \
        chk <- readr::read_tsv('~{king_file}', n_max=10); \
        if (nrow(chk) > 0) { \
            exponent <- c(5,7,9,11,13)[~{degree}]; \
            thresh <- 2^(-exponent/2); \
            kinobj <- kingToMatrix('~{king_file}', estimator='~{estimator}', thresh=thresh); \
            if (all(kinobj[upper.tri(kinobj)] > thresh)) { \
                rels <- rownames(kinobj)[2:nrow(kinobj)]; \
            } else { \
                part <- pcairPartition(kinobj, kin.thresh=thresh); \
                rels <- part[['rels']]; \
            }; \
            readr::write_tsv(tibble::tibble(FID=rels, IID=rels), 'related_samples.txt', col_names=FALSE); \
            writeLines('true', 'has_relatives.txt'); \
        } else { \
            message('No related samples found'); \
            writeLines(' ', 'related_samples.txt'); \
            writeLines('false', 'has_relatives.txt'); \
        }; \
        "
    >>>

    output {
        File related_samples = "related_samples.txt"
        Boolean has_relatives = read_boolean("has_relatives.txt")
    }

    runtime {
        docker: "uwgac/topmed-master:2.12.1"
        memory: mem_gb + " GB"
    }
}


# works on one fileset; scatter over chromosomes to apply to per-chromosome arrays
task removeSamples {
    input {
        File pgen
        File pvar
        File psam
        File samples_to_remove
        String suffix = "subset"
        Int mem_gb = 16
    }

    Int disk_size = ceil(1.5*(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB"))) + 10
    String basename = basename(pgen, ".pgen")

    command <<<
        set -e -o pipefail

        # match on IID only, so it works whether or not the psam has an FID column
        echo "#IID" > remove_iids.txt
        awk 'NF >= 2 {print $2}' ~{samples_to_remove} >> remove_iids.txt

        plink2 --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} \
            --remove remove_iids.txt \
            --output-chr chrM \
            --set-all-var-ids @:#:\$r:\$a \
            --make-pgen \
            --out ~{basename}_~{suffix}
    >>>

    output {
        File out_pgen = "~{basename}_~{suffix}.pgen"
        File out_pvar = "~{basename}_~{suffix}.pvar"
        File out_psam = "~{basename}_~{suffix}.psam"
    }

    runtime {
        docker: "quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_gb + " GB"
    }
}
