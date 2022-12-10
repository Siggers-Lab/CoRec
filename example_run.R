devtools::load_all()

corec_motifs <- run_full_analysis(
    output_directory = "./example_output",
    fluorescence_file =
        "./example_data/data_matrices/hTF_v1_SUDHL4_14jan21_br.dat",
    pbm_conditions_file = "./example_data/v1_a11_run1_exp.txt",
    # Need only one of pbm_conditions or pbm_conditions_file
    pbm_conditions = NA,
    annotation_file = "./example_data/hTF_v01_PBM_ANNOT.txt",
    reference_motifs_file =
        "./example_data/Homo_sapiens_JASPAR2022_CORE_filtered.meme",
    array_id = "v1_a11_run1",
    # I used a seed_zscore_threshold of 5 here so the example would run faster,
    #   but by default it uses 1, which is probably better
    motif_strength_threshold = 5,
    rolling_ic_threshold = 1.5,
    comparison_method = "ed",
    cluster_assignments_file = "./example_data/motif_clusters.tsv",
    pvalue_threshold = 0.05
)
