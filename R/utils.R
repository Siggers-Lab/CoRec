#' Make a data frame summarizing a list of CoRecMotifs
#'
#' Creates a data frame representation of a list of
#' [CoRecMotifs][CoRecMotif-class].
#'
#' @param corecmotifs `list`. The [CoRecMotifs][CoRecMotif-class] to summarize.
#' @param by_cluster `logical(1)`. Should the [CoRecMotifs][CoRecMotif-class] be
#'   grouped by cluster? (Default: FALSE)
#'
#' @return A data frame with information about a list of
#'   [CoRecMotifs][CoRecMotif-class].
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN")
summarize_corecmotifs <- function(corecmotifs, by_cluster = FALSE) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        assertthat::is.flag(by_cluster)
    )

    # Make sure corecmotifs is a valid list of CoRecMotifs
    corecmotifs <- check_corecmotif_list(corecmotifs)

    # Convert each CoRecMotif into a data frame
    corecmotif_df <-
        lapply(corecmotifs, as.data.frame) %>%

        # Combine all the data frames
        dplyr::bind_rows() %>%

        # Add an index mapping the data frame row to the position in the list
        dplyr::mutate(list_index = dplyr::row_number()) %>%
        dplyr::relocate(list_index, .after = array_id) %>%

        # Group replicates together
        dplyr::group_by(probe_set, pbm_condition) %>%

        # Sort by match p-value
        dplyr::arrange(match_pvalue, .by_group = TRUE) %>%

        # Add a column for the match cluster with the best p-value
        dplyr::mutate(best_match_cluster = dplyr::first(match_cluster)) %>%

        # Remove the grouping
        dplyr::ungroup()

    # Group all the CoRecMotifs from a condition that matched the same cluster
    if (by_cluster) {
        corecmotif_df <-
            corecmotif_df %>%

            # Group by match cluster and PBM condition
            dplyr::group_by(best_match_cluster, pbm_condition) %>%

            # Summarize the "best" value from each group
            dplyr::summarize(
                probe_sets =
                    paste(
                        paste(array_id, probe_set, sep = "_"),
                        collapse = ";"
                    ),
                max_motif_strength = max(motif_strength),
                max_rolling_ic = max(rolling_ic),
                match_motifs = paste(match_motif, collapse = ";"),
                min_match_pvalue = min(match_pvalue)
            )
    }

    # Return the data frame of CoRecMotif information grouped by cluster
    return(corecmotif_df)
}

#' Update the match cluster of a CoRecMotif
#'
#' Updates the `match_cluster` slot of a [CoRecMotif][CoRecMotif-class] based on
#' the provided cluster assignments and the name of the motif in the
#' CoRecMotif's `match_motif` slot.
#'
#' @param corecmotif [CoRecMotif][CoRecMotif-class]. The
#'   [CoRecMotif][CoRecMotif-class] to update.
#' @param cluster_assignments `data.frame` or `NULL`. A table mapping the
#'   reference motifs to motif clusters or NULL to reset the `match_cluster`
#'   slot to NA. See [motif_clusters] for expected columns. (Default: NULL)
#'
#' @return A [CoRecMotif][CoRecMotif-class] with its `match_cluster` slot
#'   updated.
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN")
update_cluster_match <- function(corecmotif, cluster_assignments = NULL) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        is.data.frame(cluster_assignments) || is.null(cluster_assignments)
    )

    # Make sure corecmotif is a valid CoRecMotif
    if (!methods::is(corecmotif, "CoRecMotif")) {
        stop(
            "corecmotif is not a CoRecMotif",
            call. = FALSE
        )
    }
    methods::validObject(corecmotif)

    # Clear the cluster match slot if there are no clusters or no motif match
    if (is.null(cluster_assignments) ||
        !methods::is(get_match_motif(corecmotif), "universalmotif")) {
        corecmotif@match_cluster <- NA_character_

        # Return the updated corecmotif
        return(corecmotif)
    }

    # Make sure cluster_assignments has the expected columns and remove extras
    cluster_assignments <-
        check_colnames(cluster_assignments, c("motif", "cluster"))

    # Clear the cluster match slot if the motif match isn't in the clusters
    if (!(get_match_altname(corecmotif) %in% cluster_assignments$motif)) {
        corecmotif@match_cluster <- NA_character_

        # Print a warning message
        warning(
            "Motif match altname not in cluster assignments table; ",
            "setting cluster match to NA"
        )

        # Return the updated corecmotif
        return(corecmotif)
    }

    # Figure out the cluster the best match motif is in
    best_cluster <-
        cluster_assignments %>%

        # Keep just the row corresponding to the best match motif
        dplyr::filter(motif == get_match_altname(corecmotif)) %>%

        # Pull out the cluster name for this motif
        dplyr::pull(cluster)

    # Update the cluster match slot
    corecmotif@match_cluster <- as.character(best_cluster)

    # Return the updated corecmotif
    return(corecmotif)
}

#' Update a CoRecMotif
#'
#' Updates all the values that are calculated automatically when a
#' [CoRecMotif][CoRecMotif-class] is created.
#'
#' @param corecmotif [CoRecMotif][CoRecMotif-class]. The
#'   [CoRecMotif][CoRecMotif-class] to update.
#' @param keep_match `logical(1)`. Should the `match_motif`, `match_pvalue`, and
#'   `match_cluster` slots be kept? If `FALSE`, they will be reset to `NA`.
#'   (Default: FALSE)
#'
#' @return An updated [CoRecMotif][CoRecMotif-class].
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN")
update_corecmotif <- function(corecmotif, keep_match = FALSE) {
    # Make a new CoRecMotif from the contents of the old one
    # This will update all the fields that are calculated automatically
    updated_corecmotif <-
        CoRecMotif(
            corecmotif@probe_set,
            corecmotif@pbm_condition,
            corecmotif@zscore_motif,
            corecmotif@array_id,
            corecmotif@seed_sequence
        )

    # Update the match information of the new CoRecMotif if necessary
    if (keep_match) {
        updated_corecmotif@match_motif <- corecmotif@match_motif
        updated_corecmotif@match_pvalue <- corecmotif@match_pvalue
        updated_corecmotif@match_cluster <- corecmotif@match_cluster
    }

    # Return the updated CoRecMotif
    return(updated_corecmotif)
}

