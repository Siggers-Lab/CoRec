# Load the hTFArrayAnalysis package
devtools::load_all()

# Load the PBM conditions for replicate 1
pbm_conditions_rep1 <-
    scan(
        "./example_data/hTF_v1_example_pbm_conditions_rep1.txt",
        what = character(),
        quiet = TRUE
    )

# Load the PBM conditions for replicate 2
pbm_conditions_rep2 <-
    scan(
        "./example_data/hTF_v1_example_pbm_conditions_rep2.txt",
        what = character(),
        quiet = TRUE
    )

# Set the common arguments
output_directory <- "./example_data/output"
annotation_file <- "./example_data/hTF_v1_example_annotation.tsv"
reference_motifs_file <-
    "./example_data/Homo_sapiens_JASPAR2022_CORE_filtered.meme"
cluster_assignments <-
    read.table(
        "./example_data/motif_clusters.tsv",
        header = TRUE, sep = "\t"
    )

# Make the corecmotifs for the replicate 1 array
corecmotifs_rep1 <- make_corecmotifs(
    fluorescence_file = "./example_data/hTF_v1_example_fluorescence_rep1.dat",
    pbm_conditions = pbm_conditions_rep1,
    annotation_file = annotation_file,
    output_directory = output_directory,
    output_base_name = "example_rep1",
    array_id = "v1_a11_run1"
)

# Make the corecmotifs for the replicate 2 array
corecmotifs_rep2 <- make_corecmotifs(
    fluorescence_file = "./example_data/hTF_v1_example_fluorescence_rep2.dat",
    pbm_conditions = pbm_conditions_rep2,
    annotation_file = annotation_file,
    output_directory = output_directory,
    output_base_name = "example_rep2",
    array_id = "v1_a21_run1"
)

matched_corecmotifs <-
    process_corecmotifs(
        corecmotifs = c(corecmotifs_rep1, corecmotifs_rep2),
        reference_motifs_file = reference_motifs_file,
        min_rolling_ic = 1,
        min_motif_strength = 1,
        min_n_replicates = 2,
        max_eucl_distance = 0.4,
        min_overlap = 5,
        cluster_assignments = cluster_assignments,
        max_match_pvalue = 0.05,
        meme_path = "/share/pkg.7/meme/5.3.3/install/bin/"
    )


