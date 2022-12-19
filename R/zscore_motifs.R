#' Extract corecmotif objects from z-score matrix
#'
#' Extracts a corecmotif object from the input z-score matrix for every possible
#' combination of the given seed names and PBM conditions.
#'
#' @rdname extract_corecmotifs
#'
#' @param zscore_table A data frame of PBM condition-wise fluorescence value
#'   z-scores, such as output by make_zscore_table().
#' @param pbm_condition,pbm_conditions A character vector of the column names of
#'   the PBM conditions for which to extract corecmotif objects.
#'
#' @return A list of corecmotif objects, one for each possible combination of
#'   the seed names present in \code{zscore_table} and the PBM conditions in
#'   \code{pbm_conditions}.
#' @export
#'
#' @examples
#' corec_motifs <-
#'     make_corec_motifs(
#'         zscore_table,
#'         pbm_conditions
#'     )
make_corec_motifs <-
    function(
        zscore_table,
        pbm_conditions
    ) {
        # Make a table of motif data
        motif_table <-
            zscore_table %>%

            # Keep only the seed probe row for each non-background probe set
            dplyr::filter(SNV_pos_offset == 0) %>%

            # Keep only the seed name and probe sequence columns
            dplyr::select(seed_names, probe_seq) %>%

            # Get all possible combinations of seed names and PBM conditions
            tidyr::expand_grid(pbm_conditions = pbm_conditions) %>%

            # Make the z-score motif for each seed/condition combination
            dplyr::mutate(
                zscore_motif = purrr::map2(
                    seed_names,
                    pbm_conditions,
                    make_zscore_motif,
                    zscore_table = zscore_table
                )
            )

        # Convert the data frame into a list of corecmotif objects
        corec_motifs <-
            purrr::pmap(
                list(
                    seed_name = motif_table$seed_names,
                    pbm_condition = motif_table$pbm_conditions,
                    zscore_motif = motif_table$zscore_motif,
                    seed_probe_sequence = motif_table$probe_seq
                ),
                corecmotif
            )

        # Return the list of corecmotif objects
        return(corec_motifs)
}


#' Extract corecmotif objects from z-score matrix
#'
#' @rdname extract_corecmotifs
#'
#' @param zscore_table A data frame of PBM condition-wise fluorescence value
#'   z-scores, such as output by fluorescence_to_zscore_table().
#' @param seed_name The name of the seed for which to extract a corecmotif
#'   object.
#' @param pbm_condition The column name of the PBM condition for which to
#'   extract a corecmotif object.
#'
#' @return A data frame representing the z-score motif for the seed name and PBM
#'   condition provided to \code{seed_name} and \code{pbm_condition}
#'   respectively, where the rows are nucleotides and the columns are positions
#'   in the motif.
#' @export
#'
#' @examples
#' zscore_motif <-
#'     make_zscore_motif(
#'         zscore_table,
#'         "MA0079.3_SP1",
#'         "v1_a11_run1_UT_SUDHL4_SMARCA4MIX"
#'     )
make_zscore_motif <- function(zscore_table, seed_name, pbm_condition) {
    # Get the z-score of the seed probe for this seed_name/pbm_condition combo
    seed_zscore <-
        zscore_table %>%

        # Keep only the seed probe row from the relevant probe set
        dplyr::filter(seed_names == seed_name & SNV_pos_offset == 0) %>%

        # Pull the z-score from the relevant PBM condition column
        dplyr::pull(pbm_condition)

    # Fill in a motif data frame with the seed and SV probe z-scores
    zscore_motif <-
        zscore_table %>%

        # Keep only the SV probe rows from the relevant probe set
        dplyr::filter(seed_names == seed_name & SNV_pos_offset > 0) %>%

        # Keep only the nucleotide, position, and z-score columns
        dplyr::select(
            SNV_pos_offset,
            SNV_nuc,
            zscore = !!as.symbol(pbm_condition)
        ) %>%

        # Fill in missing SNV_pos_offset/SNV_nuc combos with the seed z-score
        tidyr::complete(
            SNV_pos_offset,
            SNV_nuc,
            fill = list(zscore = seed_zscore)
        ) %>%

        # Reformat so the columns are positions and the rows are nucleotides
        tidyr::pivot_wider(
            names_from = SNV_pos_offset,
            values_from = zscore
        ) %>%

        # Make the nucleotides the row names
        tibble::column_to_rownames("SNV_nuc")

    # Return the z-score motif
    return(zscore_motif)
}

