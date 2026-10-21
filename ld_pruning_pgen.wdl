version 1.0

import "variant_filtering_pgen.wdl" as variant_tasks

workflow LD_pruning_pgen {
    input {
        # one file per chromosome; the three arrays must be the same length and in the same order
        Array[File] pgen
        Array[File] pvar
        Array[File] psam
        File? variant_file
        Int? genome_build
        Float? min_maf
        Boolean? snps_only
        Int? window_size
        Int? shift_size
        Float? r2_threshold
    }

    scatter (files in zip(zip(pgen, pvar), psam)) {
        call variant_tasks.subsetVariants {
             input:
                pgen = files.left.left,
                pvar = files.left.right,
                psam = files.right,
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
    }

    output {
        Array[File] pruned_pgen = pruneVars.out_pgen
        Array[File] pruned_pvar = pruneVars.out_pvar
        Array[File] pruned_psam = pruneVars.out_psam
    }

    meta {
        author: "Stephanie Gogarten"
        email: "sdmorris@uw.edu"
    }
}
