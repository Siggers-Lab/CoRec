# Load the CoRec package
devtools::load_all()

# Set the common arguments
reference_motifs_file <-
    "./inst/extdata/Homo_sapiens_JASPAR2022_CORE_filtered.meme"

# Load the PBM conditions for replicate 1
pbm_conditions_rep1 <-
    scan(
        "./inst/extdata/hTF_v1_example_pbm_conditions_rep1.txt",
        what = character(),
        quiet = TRUE
    )

# Load the fluorescence data for replicate 1
fluorescence_table_rep1 <-
    load_fluorescence_data(
        "./inst/extdata/hTF_v1_example_fluorescence_rep1.dat",
        pbm_conditions_rep1
    )

# Load the PBM conditions for replicate 2
pbm_conditions_rep2 <-
    scan(
        "./inst/extdata/hTF_v1_example_pbm_conditions_rep2.txt",
        what = character(),
        quiet = TRUE
    )

# Load the fluorescence data for replicate 2
fluorescence_table_rep2 <-
    load_fluorescence_data(
        "./inst/extdata/hTF_v1_example_fluorescence_rep2.dat",
        pbm_conditions_rep2
    )

# Make the CoRecMotifs for the replicate 1 array
corecmotifs_rep1 <- make_corecmotifs(
    fluorescence_table_rep1,
    fluorescence_columns = pbm_conditions_rep1,
    annotation = hTF_v1_annotation,
    array_id = "v1_a11_run1"
)

# Make the CoRecMotifs for the replicate 2 array
corecmotifs_rep2 <- make_corecmotifs(
    fluorescence_table_rep2,
    fluorescence_columns = pbm_conditions_rep2,
    annotation = hTF_v1_annotation,
    array_id = "v1_a21_run1"
)

# Filter the CoRecMotifs and match to reference motifs
matched_corecmotifs <-
    process_corecmotifs(
        corecmotifs = c(corecmotifs_rep1, corecmotifs_rep2),
        reference_motifs_file = reference_motifs_file,
        rolling_ic = 1,
        motif_strength = 1,
        n_replicates = 2,
        eucl_distance = 0.4,
        min_overlap = 5,
        cluster_assignments = motif_clusters,
        match_pvalue = 0.05,
        meme_path = "/share/pkg.7/meme/5.3.3/install/bin/"
    )

