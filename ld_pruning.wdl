version 1.0

import "variant_filtering.wdl" as variant_tasks
import "https://raw.githubusercontent.com/UW-GAC/primed-file-conversion/main/plink2_pgen2vcf.wdl" as pgen_conversion

workflow LD_pruning {
    input {
        Array[File] vcf
        File? variant_file
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
                variant_file = variant_file,
                min_maf = min_maf,
                snps_only = snps_only
        }

        call variant_tasks.pruneVars {
             input:
                pgen = subsetVariants.subset_pgen,
                pvar = subsetVariants.subset_pvar,
                psam = subsetVariants.subset_psam,
                window_size = window_size,
                shift_size = shift_size,
                r2_threshold = r2_threshold
        }

        call pgen_conversion.pgen2vcf {
            input:
                pgen = pruneVars.out_pgen,
                pvar = pruneVars.out_pvar,
                psam = pruneVars.out_psam
        }
    }

    output {
        Array[File] pruned_vcf = pgen2vcf.out_file
    }

    meta {
        author: "Stephanie Gogarten"
        email: "sdmorris@uw.edu"
    }
}
