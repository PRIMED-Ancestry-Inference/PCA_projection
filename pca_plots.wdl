version 1.0

workflow pca_plots {
    input{
        File data_file
        File? groups_file
        Int? n_pairs
    }

    call run_pca_plots {
        input: data_file = data_file, 
               groups_file = groups_file, 
               n_pairs = n_pairs
    }

    output{
        File pca_plots_pc12 = run_pca_plots.pca_plots_pc12
        Array[File] pca_plots_pairs = run_pca_plots.pca_plots_pairs
        File pca_plots_parcoord = run_pca_plots.pca_plots_parcoord
        File pca_plots = run_pca_plots.pca_plots
    }
}

task run_pca_plots {
    input{
        File data_file
        File? groups_file
        Int n_pairs = 10
    }

    command <<<
    if [ -f ~{groups_file} ]
    then
        Rscript /usr/local/PCA_projection/pca_plots.R \
            --data_file ~{data_file} \
            --groups_file ~{groups_file} \
            --n_pairs ~{n_pairs} \
            --path_to_rmd /usr/local/PCA_projection/
    else 
        Rscript /usr/local/PCA_projection/pca_plots.R \
            --data_file ~{data_file} \
            --n_pairs ~{n_pairs} \
            --path_to_rmd /usr/local/PCA_projection/
    fi
    >>>

    output{
        File pca_plots_pc12 = "out_file_pc12.png"
        Array[File] pca_plots_pairs = glob("out_file_pairs_*.png")
        File pca_plots_parcoord = "out_file_parcoord.png"
        File pca_plots = "pca_plots.html" 
    }

    runtime{
        docker: "uwgac/pca_projection:0.1.0"
    }
}