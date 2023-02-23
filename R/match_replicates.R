#' Check reproducibility of a list of CoRecMotifs
#'
#' Filter a list of [CoRecMotifs][CoRecMotif-class] based on whether they are
#' reproducible, i.e., found in multiple replicate experiments.
#'
#' @param corecmotifs the list of [CoRecMotifs][CoRecMotif-class] to filter.
#' @param n_replicates a single positive integer specifying the minimum number
#'   of replicates to require.
#' @param eucl_distance a single number specifying the maximum allowable
#'   Euclidean distance between replicate motifs or NULL to skip the replicate
#'   comparison step.
#' @param output_file the path to the RDS file where the list of reproducible
#'   [CoRecMotifs][CoRecMotif-class] will be written. If NULL (the default), no
#'   file is written.
#'
#' @return A list of [CoRecMotifs][CoRecMotif-class] that are reproducible.
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN")
match_replicates <-
    function(
        corecmotifs,
        n_replicates = 2,
        eucl_distance = 0.4,
        output_file = NULL
    ) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        assertthat::is.count(n_replicates),
        assertthat::is.number(eucl_distance) || is.null(eucl_distance),
        assertthat::is.string(output_file) || is.null(output_file)
    )

    # Make sure corecmotifs is a list
    if (!is.list(corecmotifs)) {
        corecmotifs <- list(corecmotifs)
    }

    # Make a dataframe summarizing the CoRecMotifs
    corecmotif_df <-
        summarize_corecmotifs(corecmotifs) %>%

        # Add a column with the row number to map to the list of CoRecMotifs
        dplyr::mutate(index = dplyr::row_number()) %>%

        # Group replicates together
        dplyr::group_by(seed_name, pbm_condition) %>%

        # Remove any replicate groups that don't have enough motifs
        dplyr::filter(dplyr::n() >= n_replicates)

    # Make a list of lists of replicate CoRecMotifs
    grouped_corecmotifs <-
        corecmotif_df %>%

        # Each internal list is all the motifs that are replicates of each other
        dplyr::group_map(~ c(corecmotifs[.x$index]), .keep = TRUE)

    # If not filtering by similarity, return the filtered list now
    if (is.null(eucl_distance)) {
        # Flatten the list before returning
        grouped_corecmotifs <-
            purrr::flatten(grouped_corecmotifs)

        # Try to save the replicated CoRecMotifs as an RDS file if necessary
        try_catch_save_output(grouped_corecmotifs, output_file, "rds")

        return(grouped_corecmotifs)
    }

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
            dplyr::filter(rowname != name & value <= eucl_distance) %>%

            # Get the names of all the remaining similar motifs
            dplyr::pull(rowname)

        # Return the replicate CoRecMotifs that are replicating
        return(group[motif_names %in% replicated_motif_names])
    })

    # Remove any groups that don't have enough motifs
    replicated_motifs <-
        lapply(replicated_motifs, function(group) {
            if (length(group) < n_replicates) {
                return(NULL)
            }
            return(group)
        }) %>%

        # Remove the grouping so it's a flat list of CoRecMotifs again
        purrr::flatten()

    # Try to save the replicated CoRecMotifs as an RDS file if necessary
    try_catch_save_output(replicated_corecmotifs, output_file, "rds")

    # Return the list of replicated CoRecMotifs
    return(replicated_motifs)
}

