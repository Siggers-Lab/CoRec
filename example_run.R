devtools::load_all()

# Load the PBM conditions from a file
pbm_conditions <-
    scan(
        "./example_data/hTF_v1_example_pbm_conditions.txt",
        what = character(),
        quiet = TRUE
    )

# Set all the arguments
output_directory <- "./example_data/example_output"
fluorescence_file <- "./example_data/hTF_v1_example_fluorescence.dat"
annotation_file <- "./example_data/hTF_v1_example_annotation.tsv"
reference_motifs_file <-
    "./example_data/Homo_sapiens_JASPAR2022_CORE_filtered.meme"
output_base_name <- "example"
array_id <- "v1_a11_run1"
motif_strength_threshold <- 1
rolling_ic_threshold <- 1.5
comparison_method <- "ed"
cluster_assignments_file <- "./example_data/motif_clusters.tsv"
pvalue_threshold <- 0.05

# Run the analysis
corec_motifs <- make_corecmotifs(
    fluorescence_file = fluorescence_file,
    pbm_conditions = pbm_conditions,
    annotation_file = annotation_file,
    output_directory = output_directory,
    output_base_name = output_base_name,
    array_id = array_id
)
