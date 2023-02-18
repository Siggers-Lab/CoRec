match_replicates <-
    function(
        corecmotifs,
        min_n_replicates = 2,
        max_eucl_distance = 0.4
    ) {
    # Make a dataframe summarizing the corecmotifs
    corecmotif_df <-
        summarize_corecmotifs(corecmotifs) %>%

        # Add a column with the row number to map to the list of corecmotifs
        dplyr::mutate(index = dplyr::row_number()) %>%

        # Group replicates together
        dplyr::group_by(seed_name, pbm_condition) %>%

        # Remove any replicate groups that don't have enough motifs
        dplyr::filter(dplyr::n() >= min_n_replicates)

    # Make a list of lists of replicate corecmotifs
    grouped_corecmotifs <-
        corecmotif_df %>%

        # Each internal list is all the motifs that are replicates of each other
        dplyr::group_map(~ c(corecmotifs[.x$index]), .keep = TRUE)

    # Make sure replicate motifs are actually replicating (i.e., similar)
    replicated_motifs <- lapply(grouped_corecmotifs, function(group) {
        # Get a list of the PPMs in this group of replicates
        ppms <- lapply(group, get_ppm)

        # Get the names of the motifs in this group
        motif_names <- lapply(group, get_motif_name)

        # Compare all the motifs to each other
        motif_comparison <-
            universalmotif::compare_motifs(
                ppms,
                method = "EUCL",
                # Don't bother with the reverse complements
                tryRC = FALSE,
                # Align all the columns
                min.overlap = 100,
                min.mean.ic = 0
            )

        # Get the names of similar motifs
        replicated_motif_names <-
            motif_comparison %>%

            # Convert the distance matrix to a dataframe
            as.data.frame() %>%

            # Convert the rownames into a column
            tibble::rownames_to_column() %>%

            # Convert to long format (columns: motif 1, motif 2, distance)
            tidyr::pivot_longer(cols = colnames(motif_comparison)) %>%

            # Filter out self comparisons and dissimilar motifs
            dplyr::filter(rowname != name & value <= max_eucl_distance) %>%

            # Get the names of all the remaining similar motifs
            dplyr::pull(rowname)

        # Return the replicate corecmotifs that are replicating
        return(group[motif_names %in% replicated_motif_names])
    })

    # Remove any groups that don't have enough motifs
    replicated_motifs <-
        lapply(replicated_motifs, function(group) {
            if (length(group) < min_n_replicates) {
                return(NULL)
            }
            return(group)
        }) %>%

        # Remove the grouping so it's a flat list of corecmotifs again
        purrr::flatten()

    # Return the list of replicated corecmotifs
    return(replicated_motifs)
}

