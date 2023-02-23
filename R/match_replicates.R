#' Check reproducibility of a list of \linkS4class{CoRecMotif} objects
#'
#' Filter a list of \linkS4class{CoRecMotif} objects based on whether they are
#' reproducible, i.e., found in multiple replicate experiments.
#'
#' @param corecmotifs the list of CoRecMotif objects to filter.
#' @param min_n_replicates a single positive integer specifying the minimum
#'   number of replicates to require.
#' @param max_eucl_distance a single number specifying the maximum allowable
#'   Euclidean distance between replicate motifs or NULL to skip the replicate
#'   comparison step.
#'
#' @return A list of \linkS4class{CoRecMotif} objects that are reproducible.
#'
#' @export
#'
#' @examples
#' # Load example CoRecMotifs
#' corecmotifs_rep1 <-
#'     readRDS(
#'         "example_data/output/example_rep1_v1_a11_run1_all_corecmotifs.rds"
#'     )
#' corecmotifs_rep2 <-
#'     readRDS(
#'         "example_data/output/example_rep2_v1_a21_run1_all_corecmotifs.rds"
#'     )
#'
#' # Filter out dissimilar "replicate" motifs
#' replicated_corecmotifs <-
#'     match_replicates(c(corecmotifs_rep1, corecmotifs_rep2))
match_replicates <-
    function(
        corecmotifs,
        min_n_replicates = 2,
        max_eucl_distance = 0.4
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
        dplyr::filter(dplyr::n() >= min_n_replicates)

    # Make a list of lists of replicate CoRecMotifs
    grouped_corecmotifs <-
        corecmotif_df %>%

        # Each internal list is all the motifs that are replicates of each other
        dplyr::group_map(~ c(corecmotifs[.x$index]), .keep = TRUE)

    # If not filtering by similarity, return the filtered list now
    if (is.null(max_eucl_distance)) {
        # Flatten the list before returning
        grouped_corecmotifs <-
            purrr::flatten(grouped_corecmotifs)

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

