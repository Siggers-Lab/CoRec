#' Compare CoRecMotifs across different PBM conditions
#'
#' Compares the [CoRecMotifs][CoRecMotif-class] in each PBM condition to the
#' CoRecMotifs from the same probe set in all the other PBM conditions.
#'
#' @inheritParams filter_corecmotifs
#' @param corecmotifs `list`. The [CoRecMotifs][CoRecMotif-class] to compare to
#'   each other.
#' @param pbm_conditions `character`. The names of the individual PBM conditions
#'   to compare.
#' @param pbm_conditions_group `character(1)` or `NULL`. The name of the group
#'   of PBM conditions to compare. This name will be combined with the probe set
#'   ID and a number to identify groups of motifs that are similar across
#'   conditions. (Default: NULL)
#' @param eucl_distance `numeric(1)`. The maximum allowable Euclidean distance
#'   between conditions to group. (Default: 0.4)
#'
#' @return A data frame with comparison information about a list of
#'   [CoRecMotifs][CoRecMotif-class].
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN")
compare_conditions <-
    function(
        corecmotifs,
        pbm_conditions,
        pbm_conditions_group = NULL,
        eucl_distance = 0.25,
        check_corecmotifs = TRUE
    ) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        is.character(pbm_conditions),
        assertthat::is.string(pbm_conditions_group) ||
            is.null(pbm_conditions_group),
        assertthat::is.number(eucl_distance),
        assertthat::is.flag(check_corecmotifs)
    )

    # Summarize the list of all CoRecMotifs
    corecmotif_df <-
        summarize_corecmotifs(
            corecmotifs,
            check_corecmotifs = check_corecmotifs
        )

    # Make a name for this group of conditions if one is not provided
    if (is.null(pbm_conditions_group) || pbm_conditions_group == "") {
        pbm_conditions_group <- paste(pbm_conditions, collapse = "_")
    }

    # Keep only the CoRecMotifs from the relevant PBM conditions
    matching_corecmotifs <-
        filter_corecmotifs(
            corecmotifs,
            pbm_condition = pbm_conditions,
            check_corecmotifs = FALSE
        )

    # Summarize just the matching CoRecMotifs
    matching_corecmotif_df <-
        summarize_corecmotifs(
            matching_corecmotifs,
            check_corecmotifs = FALSE
        ) %>%

        # Group all the motifs from the same probe set together
        dplyr::group_by(probe_set)

    # Make a list of lists of CoRecMotifs from the same probe set
    grouped_corecmotifs <-
        matching_corecmotif_df %>%

        # Each internal list is all the motifs from the same probe set
        dplyr::group_map(~ c(matching_corecmotifs[.x$list_index]), .keep = TRUE)

    # Compare motifs from the same probe set in different PBM conditions
    motif_comparisons <- lapply(grouped_corecmotifs, function(group) {
        # Get the name of the probe set
        probe_set <- get_probe_set(group[[1]])

        # Get a list of all the PBM conditions for this probe set
        group_pbm_conditions <-
            vapply(group, get_pbm_condition, character(1)) %>%

            unique()

        # If there's only one PBM condition, there's nothing to compare
        if (length(group_pbm_conditions) == 1) {
            motif_comparison <-
                data.frame(
                    "probe_set" = probe_set,
                    "pbm_condition" = group_pbm_conditions,
                    "group" = paste0(pbm_conditions_group, "_", probe_set, "_1")
                )
            return(motif_comparison)
        }

        # Get a list of the PPMs from this probe set
        motifs <- lapply(group, get_motif)

        # Get the names of the motifs in this group
        motif_names <- vapply(group, get_motif_name, character(1))

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

        # Find the distance between each condition
        motif_comparison <-
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

            # Add the information about motif_1
            dplyr::left_join(
                corecmotif_df,
                by = c("motif_name_1" = "motif_name")
            ) %>%

            # Add the information about motif_2
            dplyr::left_join(
                corecmotif_df,
                by = c("motif_name_2" = "motif_name", "probe_set"),
                suffix = c("_1", "_2")
            ) %>%

            # Group by the PBM conditions being compared
            dplyr::group_by(probe_set, pbm_condition_1, pbm_condition_2) %>%

            # Take the minimum distance between two conditions
            dplyr::summarise(min_distance = min(distance))

        # Group similar conditions together
        group_assignments <-
            motif_comparison %>%

            # Convert to wide format with between condition distances
            tidyr::pivot_wider(
                id_cols = "pbm_condition_1",
                names_from = "pbm_condition_2",
                values_from = "min_distance"
            ) %>%

            # Convert the PBM condition column into the row names
            tibble::column_to_rownames("pbm_condition_1") %>%

            # Convert to a dist object
            stats::as.dist() %>%

            # Cluster the conditions using single linkage
            stats::hclust(method = "single") %>%

            # Cut into groups that are separated by at least eucl_distance
            stats::cutree(h = eucl_distance)

        # Make a data frame of the group assignments
        group_assignments <-
            data.frame(
                "probe_set" = rep(probe_set, length(group_assignments)),
                "pbm_condition" = names(group_assignments),
                "group" =
                    paste(
                        pbm_conditions_group,
                        probe_set,
                        group_assignments,
                        sep = "_"
                    )
            )

        # Return the data frame of group assignments
        return(group_assignments)
    })

    motif_comparison_df <-
        motif_comparisons %>%

        # Combine all the data frames
        dplyr::bind_rows() %>%

        # Add the motif information
        dplyr::left_join(
            corecmotif_df,
            by = c("probe_set", "pbm_condition")
        ) %>%

        # Group similar PBM conditions for the same probe set together
        dplyr::group_by(probe_set, group) %>%

        # Sort by match p-value
        dplyr::arrange(match_pvalue, .by_group = TRUE) %>%

        # Add a column for the match cluster with the best p-value
        dplyr::mutate(group_match_cluster = dplyr::first(match_cluster)) %>%

        # Remove the grouping
        dplyr::ungroup()
}

