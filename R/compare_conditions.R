#' Compare CoRecMotifs across different PBM conditions
#'
#' Compares the [CoRecMotifs][CoRecMotif-class] in each PBM condition to the
#' CoRecMotifs from the same probe set in all the other PBM conditions.
#'
#' @param corecmotifs `list`. The [CoRecMotifs][CoRecMotif-class] to compare to
#'   each other.
#'
#' @return A data frame with comparison information about a list of
#'   [CoRecMotifs][CoRecMotif-class].
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN")
compare_conditions <- function(corecmotifs) {
    # Make a data frame summarizing the CoRecMotifs
    corecmotif_df <-
        # summarize_corecmotifs(matching_corecmotifs) %>%
        summarize_corecmotifs(corecmotifs) %>%

        # Group all the motifs from the same probe set together
        dplyr::group_by(probe_set)

    # Make a list of lists of CoRecMotifs from the same probe set
    grouped_corecmotifs <-
        corecmotif_df %>%

        # Each internal list is all the motifs from the same probe set
        dplyr::group_map(~ c(corecmotifs[.x$list_index]), .keep = TRUE)

    # Compare motifs from the same probe set in different PBM conditions
    motif_comparisons <- lapply(grouped_corecmotifs, function(group) {
        # Get a list of the PPMs in this group of replicates
        motifs <- lapply(group, get_motif)

        # Get the names of the motifs in this group
        motif_names <- vapply(group, get_motif_name, character(1))

        # If there's only one motif in this group, there's nothing to compare
        if (length(motif_names) == 1) {
            motif_comparison <-
                data.frame(
                    "motif_name_1" = motif_names,
                    "motif_name_2" = motif_names,
                    "distance" = 0
                )
            return(motif_comparison)
        }

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

        # Convert the distance matrix to a data frame
        motif_comparison <-
            motif_comparison %>%

            as.data.frame() %>%

            # Convert the rownames into a column
            tibble::rownames_to_column("motif_name_1") %>%

            # Convert to long format
            tidyr::pivot_longer(
                cols = colnames(motif_comparison)[-1],
                names_to = "motif_name_2",
                values_to = "distance"
            )

        # Return the data frame of between condition distances
        return(motif_comparison)
    })

    motif_comparison_df <-
        motif_comparisons %>%

        # Combine all the data frames
        dplyr::bind_rows() %>%

        # Add the information about motif_1
        dplyr::left_join(
            corecmotif_df,
            by = c("motif_name_1" = "motif_name")
        ) %>%

        # Add the information about motif_2
        dplyr::left_join(
            corecmotif_df,
            by = c("motif_name_2" = "motif_name", "probe_set", "seed_sequence"),
            suffix = c("_1", "_2")
        ) %>%

        dplyr::filter(pbm_condition_1 <= pbm_condition_2) %>%

        # Group by the PBM conditions being compared
        dplyr::group_by(probe_set, pbm_condition_1, pbm_condition_2) %>%

        # Take the average distance
        dplyr::summarise(
            mean_distance = mean(distance),
            motif_names_1 = paste(unique(motif_name_1), collapse = ";"),
            motif_names_2 = paste(unique(motif_name_2), collapse = ";"),
            list_indices_1 = paste(unique(list_index_1), collapse = ";"),
            list_indices_2 = paste(unique(list_index_2), collapse = ";"),
            best_match_cluster_1 = unique(best_match_cluster_1),
            best_match_cluster_2 = unique(best_match_cluster_2),
            best_match_pvalue_1 = ifelse(
                all(is.na(match_pvalue_1)),
                NA,
                min(stats::na.omit(match_pvalue_1))
            ),
            best_match_pvalue_2 = ifelse(
                all(is.na(match_pvalue_2)),
                NA,
                min(stats::na.omit(match_pvalue_2))
            )
        )
}
