version 1.0

workflow pca_plots {
    input{
        File data_file
        File groups_file
        Int n_pairs
    }

    call run_pca_plots {
        input: data_file = data_file, 
               groups_file = groups_file, 
               n_pairs = n_pairs
    }

    output{
        File pca_plots_pc12 = run_pca_plots.pca_plots_pc12
        File pca_plots_pairs = run_pca_plots.pca_plots_pairs
        File pca_plots_parcoord = run_pca_plots.pca_plots_parcoord
    }
}

task run_pca_plots {
    input{
        File data_file
        File groups_file
        Int n_pairs
    }

    command <<<
    Rscript /usr/local/PCA_projection/pca_plots.R \
    --data_file ~{data_file} \ 
    --groups_file ~{groups_file} \
    --n_pairs ~{n_pairs}
    >>>

    output{
        File pca_plots_pc12 = "out_file_pc12.pdf"
        File pca_plots_pairs = "out_file_pairs.pdf"
        File pca_plots_parcoord = "out_file_parcoord.pdf"
    }

    runtime{
        docker: "primed_ancestry_inference/pca_projection:0.5.0"
    }
}