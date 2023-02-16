#' Create \linkS4class{corecmotif} objects from a table of z-score data
#'
#' Creates a list of \linkS4class{corecmotif} objects for all possible
#' combinations of the probe sets present in \code{zscore_table} and the PBM
#' conditions given in \code{zscore_columns}.
#'
#' @param zscore_table a data frame of z-scores and annotations for each probe.
#'   See 'Details' of \code{\link{make_fluorescence_table}} for a description of
#'   the expected annotation columns.
#' @param zscore_columns a character vector specifying the names of the columns
#'   of \code{zscore_table} that contain z-score data.
#' @param output_file the name of the RDS file where the list of corecmotif
#'   objects will be written. If NULL (the default), no file is written.
#'
#' @return A list of \linkS4class{corecmotif} objects, one for each possible
#'   combination of the probe sets in \code{zscore_table} and the PBM conditions
#'   listed in \code{zscore_columns}.
#'
#' @export
#'
#' @examples
#' # Load the example z-score table
#' zscore_table <-
#'     read.table(
#'         "example_data/example_output/example_v1_a11_run1_zscores.tsv",
#'         header = TRUE,
#'         sep = "\t",
#'         stringsAsFactors = FALSE
#'     )
#'
#' # Make a list of corecmotif objects
#' corec_motifs <-
#'     make_corec_motifs(
#'         zscore_table,
#'         zscore_columns = c(
#'             "v1_a11_run1_UT_SUDHL4_SMARCA4MIX",
#'             "v1_a11_run1_UT_SUDHL4_HDAC1MIX",
#'             "v1_a11_run1_UT_SUDHL4_SUZ12",
#'             "v1_a11_run1_UT_SUDHL4_PRMT5"
#'         )
#'     )
make_corec_motifs <-
    function(
        zscore_table,
        zscore_columns,
        output_file = NULL
    ) {
    # Make a table of motif data
    motif_table <-
        zscore_table %>%

        # Keep only the seed probe row for each non-background probe set
        dplyr::filter(SNV_pos_offset == 0) %>%

        # Keep only the probe set name and probe sequence columns
        dplyr::select(seed_names, probe_seq) %>%

        # Get all possible combinations of probe set names and PBM conditions
        tidyr::expand_grid(pbm_conditions = zscore_columns) %>%

        # Make the z-score motif for each probe set/condition combination
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

    # Save the list of all corecmotifs as an RDS file if necessary
    if (!is.null(output_file)) {
        saveRDS(
            corec_motifs,
            output_file
        )
    }

    # Return the list of corecmotif objects
    return(corec_motifs)
}

#' Create a z-score motif for a given probe set and PBM condition
#'
#' Creates a z-score motif for a given probe set and PBM condition. The rows
#' correspond to nucleotides and the columns correspond to positions in the
#' motif.
#'
#' @param zscore_table a data frame of z-scores and annotations for each probe.
#'   See 'Details' of \code{\link{make_fluorescence_table}} for a description of
#'   the expected annotation columns.
#' @param probe_set a character string containing the name of the probe set
#'   for which to create the z-score motif.
#' @param pbm_condition a character string containing the name of the PBM
#'   condition for which to create the z-score motif.
#'
#' @return A data frame with rows corresponding to nucleotides and columns
#'   corresponding to positions in the motif. Each cell is filled with the
#'   z-score of the relevant probe in the given PBM condition.
#'
#' @export
#'
#' @examples
#' # Load the example z-score table
#' zscore_table <-
#'     read.table(
#'         "example_data/example_output/example_v1_a11_run1_zscores.tsv",
#'         header = TRUE,
#'         sep = "\t",
#'         stringsAsFactors = FALSE
#'     )
#'
#' # Create a z-score motif
#' zscore_motif <-
#'     make_zscore_motif(
#'         zscore_table,
#'         probe_set = "MA0052.3_MEF2A",
#'         pbm_condition = "v1_a11_run1_UT_SUDHL4_PRMT5"
#'     )
make_zscore_motif <- function(zscore_table, probe_set, pbm_condition) {
    # Get the z-score of the seed probe for this seed_name/pbm_condition combo
    seed_zscore <-
        zscore_table %>%

        # Keep only the seed probe row from the relevant probe set
        dplyr::filter(seed_names == probe_set & SNV_pos_offset == 0) %>%

        # Pull the z-score from the relevant PBM condition column
        dplyr::pull(pbm_condition)

    # Fill in a motif data frame with the seed and SV probe z-scores
    zscore_motif <-
        zscore_table %>%

        # Keep only the SV probe rows from the relevant probe set
        dplyr::filter(seed_names == probe_set & SNV_pos_offset > 0) %>%

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

