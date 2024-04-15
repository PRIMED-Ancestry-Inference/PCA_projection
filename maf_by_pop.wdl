version 1.0

workflow select_variants_by_pop_maf {
    input {
        Array[File] vcf
        Float min_maf
        String population_descriptor
        Array[String]? population_labels
        String workspace_name
        String workspace_namespace
        Boolean? snps_only
    }

    if (!defined(population_labels)) {
        call find_pop_labels {
            input:
                population_descriptor = population_descriptor,
                workspace_name = workspace_name,
                workspace_namespace = workspace_namespace
        }
    }

    scatter (pop in select_first([population_labels, find_pop_labels.labels])) {
        call samples_in_pop {
            input:
                population_descriptor = population_descriptor,
                population_label = pop,
                workspace_name = workspace_name,
                workspace_namespace = workspace_namespace
        }

        scatter (file in vcf) {
            call maf_by_pop {
                input:
                    vcf = file,
                    samples = samples_in_pop.samples,
                    min_maf = min_maf,
                    snps_only = snps_only
            }
        }
    }

    call combine_variants {
        input:
            variants = flatten(maf_by_pop.selected_variants)
    }

    output {
        File maf_filtered_variants = combine_variants.combined_variants
    }
}


task find_pop_labels {
    input {
        String population_descriptor
        String workspace_name
        String workspace_namespace
    }

    command <<<
        Rscript -e "\
        library(AnVIL); \
        library(dplyr); \
        library(tidyr); \
        pop <- avtable('population_descriptor', name='~{workspace_name}', namespace='~{workspace_namespace}'); \
        dat <- separate_longer_delim(pop, starts_with('population'), delim='|'); \
        dat <- mutate(dat, across(starts_with('population'), stringr::str_trim)); \
        dat <- pivot_wider(dat, names_from=population_descriptor, values_from=population_label); \
        stopifnot(is.element('~{population_descriptor}', names(dat))); \
        labels <- unlist(distinct(dat, ~{population_descriptor})); \
        writeLines(labels, 'labels.txt'); \
        "
    >>>

    output {
        Array[String] labels = read_lines("labels.txt")
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-gcr-public/anvil-rstudio-bioconductor:3.17.0"
    }
}


task samples_in_pop {
    input {
        String population_descriptor
        String population_label
        String workspace_name
        String workspace_namespace
    }

    command <<<
        Rscript -e "\
        library(AnVIL); \
        library(dplyr); \
        library(tidyr); \
        pop <- avtable('population_descriptor', name='~{workspace_name}', namespace='~{workspace_namespace}'); \
        samp <- avtable('sample', name='~{workspace_name}', namespace='~{workspace_namespace}'); \
        dat <- select(samp, 'subject_id', 'sample_id'); \
        dat <- inner_join(dat, pop); \
        dat <- separate_longer_delim(dat, starts_with('population'), delim='|'); \
        dat <- mutate(dat, across(starts_with('population'), stringr::str_trim)); \
        dat <- pivot_wider(dat, names_from=population_descriptor, values_from=population_label); \
        stopifnot(is.element('~{population_descriptor}', names(dat))); \
        dat <- filter(dat, is.element(~{population_descriptor}, '~{population_label}') ); \
        samples <- dat[['sample_id']]; \
        writeLines(samples, 'samples.txt'); \
        "
    >>>

    output {
        File samples = "samples.txt"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-gcr-public/anvil-rstudio-bioconductor:3.17.0"
    }
}


task maf_by_pop {
    input {
        File vcf
        File samples
        Float min_maf
        Boolean snps_only = false
        Int mem_gb = 8
    }

    Int disk_size = ceil(2.5*(size(vcf, "GB"))) + 5
    String filename = basename(vcf)
    String prefix = if (sub(filename, ".bcf", "") != filename) then "--bcf" else "--vcf"

    command <<<
        /plink2 ~{prefix} ~{vcf} \
            --keep ~{samples} \
            --maf ~{min_maf} \
            --set-missing-var-ids @:#:\$r:\$a \
            --rm-dup exclude-all \
            ~{true="--snps-only 'just-acgt'" false="" snps_only} \
            --write-snplist --out maf_filter
    >>>

    output {
        File selected_variants = "maf_filter.snplist"
    }

    runtime {
        docker: "quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_gb + " GB"
    }
}


task combine_variants {
    input {
        Array[File] variants
    }

    Int disk_size = ceil(1.5*(size(variants, "GB"))) + 5

    command <<<
        Rscript -e "\
        files <- readLines('~{write_lines(variants)}'); \
        vars <- unique(unlist(lapply(files, readLines))); \
        writeLines(vars, 'variants.txt')
        "
    >>>

    output {
        File combined_variants = "variants.txt"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-gcr-public/anvil-rstudio-bioconductor:3.17.0"
        disks: "local-disk " + disk_size + " SSD"
        memory: disk_size + " GB"
    }
}
