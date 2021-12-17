devtools::load_all()

corec_motifs <- run_full_analysis(
    output_directory = "./example_output",
    matrix_directory = "./example_data/data_matrices",
    matrix_base_name = "hTF_v1_SUDHL4_14jan21",
    pbm_conditions_file = "./example_data/v1_a11_run1_exp.txt",
    pbm_conditions = NA, # Need only one of pbm_conditions/pbm_conditions_file
    annotation_file = "./example_data/hTF_v01_PBM_ANNOT.txt",
    reference_motifs_file = "./example_data/JASPAR2018_hTF_only.meme",
    run_tag = "v1_a11_run1",
    score_method = "seed_zscore",
    score_threshold = 0.5,
    comparison_method = "ed",
    pvalue_threshold = 0.01
)
