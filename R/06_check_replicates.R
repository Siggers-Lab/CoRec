#' Check reproducibility of a list of CoRecMotifs
#'
#' Filter a list of [CoRecMotifs][CoRecMotif-class] based on whether they are
#' reproducible, i.e., found in multiple replicate experiments.
#'
#' This function is intended to be used in concert with [filter_corecmotifs()].
#' For example, you may want to remove any PBM condition/probe set combination
#' that does not have at least 3 replicate [CoRecMotifs][CoRecMotif-class] with
#' a motif strength of at least 1. In this case you would use
#' [filter_corecmotifs()] with `motif_strength = 1` to remove individual
#' [CoRecMotifs][CoRecMotif-class] with low motif strengths followed by
#' `check_replicates()` with `n_replicates = 3` to remove any replicate groups
#' that have fewer than 3 remaining [CoRecMotifs][CoRecMotif-class].
#'
#' @inheritParams annotate_fluorescence_table
#' @inheritParams filter_corecmotifs
#' @param n_replicates `integer(1)`. The minimum number of replicates to
#'   require. Set this to 1 to skip filtering by number of replicates. (Default:
#'   2)
#' @param eucl_distance `numeric(1)` or `NULL`. The maximum allowable Euclidean
#'   distance between replicate motifs or NULL to skip the replicate comparison
#'   step. (Default: 0.4)
#'
#' @return A list of [CoRecMotifs][CoRecMotif-class] that are replicating.
#'
#' @seealso [filter_corecmotifs()] for filtering individual
#'   [CoRecMotifs][CoRecMotif-class].
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN")
check_replicates <-
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

    # Make a data frame summarizing the CoRecMotifs
    corecmotif_df <-
        summarize_corecmotifs(corecmotifs) %>%

        # Group replicates together
        dplyr::group_by(probe_set, pbm_condition) %>%

        # Remove any replicate groups that don't have enough motifs
        dplyr::filter(dplyr::n() >= n_replicates)

    # Make a list of lists of replicate CoRecMotifs
    grouped_corecmotifs <-
        corecmotif_df %>%

        # Each internal list is all the motifs that are replicates of each other
        dplyr::group_map(~ c(corecmotifs[.x$list_index]), .keep = TRUE)

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
        # If there's only one motif in this group, there's nothing to compare
        if (length(group) == 1) {
            return(group)
        }

        # Get a list of the PPMs in this group of replicates
        motifs <- lapply(group, get_motif)

        # Get the names of the motifs in this group
        motif_names <- lapply(group, get_motif_name)

        # Make sure there aren't any duplicate names
        if (any(duplicated(motif_names))) {
            warning(
                "CoRecMotif names are not unique!\n",
                "This could cause unexpected behavior. Please make sure the ",
                "names are unique by providing a different array ID for each ",
                "motif in a group of replicates.",
                call. = FALSE
            )
        }

        # Compare all the motifs to each other
        motif_comparison <-
            universalmotif::compare_motifs(
                motifs,
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

            # Convert the distance matrix to a data frame
            as.data.frame() %>%

            # Convert the row names into a column
            tibble::rownames_to_column("motif_name_1") %>%

            # Convert to long format
            tidyr::pivot_longer(
                cols = colnames(motif_comparison),
                names_to = "motif_name_2",
                values_to = "distance"
            ) %>%

            # Filter out self comparisons and dissimilar motifs
            dplyr::filter(
                motif_name_1 != motif_name_2 & distance <= eucl_distance
            ) %>%

            # Get the names of all the remaining similar motifs
            dplyr::pull(motif_name_1)

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
    try_catch_save_output(replicated_motifs, output_file, "rds")

    # Return the list of replicated CoRecMotifs
    return(replicated_motifs)
}

