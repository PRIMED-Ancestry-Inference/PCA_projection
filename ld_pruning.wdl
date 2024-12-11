version 1.0

import "variant_filtering.wdl" as variant_tasks
import "https://raw.githubusercontent.com/UW-GAC/primed-file-conversion/main/plink2_bed2vcf.wdl" as bed_conversion

workflow LD_pruning {
    input {
        Array[File] vcf
        File? variant_file
        Int? genome_build
        Float? min_maf
        Boolean? snps_only
        Int? window_size
        Int? shift_size
        Float? r2_threshold
    }

    scatter (file in vcf) {
        call variant_tasks.subsetVariants {
             input:
                vcf = file,
                variant_files = select_all([variant_file]),
                genome_build = genome_build,
                min_maf = min_maf,
                snps_only = snps_only
        }

        call variant_tasks.pruneVars {
             input:
                bed = subsetVariants.subset_bed,
                bim = subsetVariants.subset_bim,
                fam = subsetVariants.subset_fam,
                window_size = window_size,
                shift_size = shift_size,
                r2_threshold = r2_threshold
        }

        call bed_conversion.bed2vcf {
            input:
                bed_file = pruneVars.out_bed,
                bim_file = pruneVars.out_bim,
                fam_file = pruneVars.out_fam
        }
    }

    output {
        Array[File] pruned_vcf = bed2vcf.out_file
    }

    meta {
        author: "Stephanie Gogarten"
        email: "sdmorris@uw.edu"
    }
}
