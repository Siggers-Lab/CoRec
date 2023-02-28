# Load the CoRec package
devtools::load_all()

# Set up some arguments
reference_motifs_file <-
    "./inst/extdata/Homo_sapiens_JASPAR2022_CORE_filtered.meme"

pbm_conditions <-
    c(
        "UT_SUDHL4_SWISNF_mix",
        "UT_SUDHL4_HDAC_mix",
        "UT_SUDHL4_PRMT5",
        "UT_SUDHL4_JMJD2A"
    )

# Use the example fluorescence data table for replicate 1
fluorescence_table_1 <- example_fluorescence_table

# Load the fluorescence data for replicate 2
fluorescence_table_2 <-
    load_fluorescence_data(
        "./inst/extdata/example_fluorescence_data_2.dat",
        pbm_conditions = pbm_conditions
    )

# Make the CoRecMotifs for the replicate 1 array
corecmotifs_1 <- make_corecmotifs(
    fluorescence_table_1,
    fluorescence_columns = pbm_conditions,
    annotation = hTF_v1_annotation,
    array_id = "v1_a11_run1"
)

# Make the CoRecMotifs for the replicate 2 array
corecmotifs_2 <- make_corecmotifs(
    fluorescence_table_2,
    fluorescence_columns = pbm_conditions,
    annotation = hTF_v1_annotation,
    array_id = "v1_a21_run1"
)

# Filter the CoRecMotifs and match to reference motifs
# This may take a while
matched_corecmotifs <-
    process_corecmotifs(
        corecmotifs = c(corecmotifs_1, corecmotifs_2),
        reference_motifs_file = reference_motifs_file,
        rolling_ic = 1,
        motif_strength = 5,
        n_replicates = 2,
        eucl_distance = 0.4,
        min_overlap = 5,
        cluster_assignments = motif_clusters,
        match_pvalue = 0.05,
        meme_path = "/share/pkg.7/meme/5.3.3/install/bin/"
    )

# Plot the delta z-score motif next to the ICM of the matching reference motif
plot_corecmotif(matched_corecmotifs[[1]])

# You can change the type of logo to plot with the *_logo_type arguments
plot_corecmotif(
    matched_corecmotifs[[1]],
    corecmotif_logo_type = "ICM",
    reference_logo_type = "none"
)

# You can generate plots for all the matched CoRecMotifs with lapply
corecmotif_plots <- lapply(matched_corecmotifs, plot_corecmotif)

# Display them all in one plot with cowplot::plot_grid()
cowplot::plot_grid(plotlist = corecmotif_plots, ncol = 2)

