#' Create CoRecMotifs from a table of z-score data
#'
#' Creates a list of [CoRecMotifs][CoRecMotif-class] for all possible
#' combinations of the probe sets present in `zscore_table` and the PBM
#' conditions given in `zscore_columns`.
#'
#' @param zscore_table a data frame of z-scores and annotations for each probe.
#'   See 'Details' of [annotate_fluorescence_table()] for a description of the
#'   expected annotation columns.
#' @param zscore_columns a character vector specifying the names of the columns
#'   of `zscore_table` that contain z-score data.
#' @param array_id an optional (but recommended) tag specifying the particular
#'   array/experiment the fluorescence data is from. (Default: NULL)
#' @param output_file the path to the RDS file where the list of
#'   [CoRecMotifs][CoRecMotif-class] will be written. If NULL, no file is
#'   written. (Default: NULL)
#'
#' @return A list of [CoRecMotifs][CoRecMotif-class], one for each possible
#'   combination of the probe sets in `zscore_table` and the PBM conditions
#'   listed in `zscore_columns`.
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN")
zscore_table_to_corecmotifs <-
    function(
        zscore_table,
        zscore_columns,
        array_id = NULL,
        output_file = NULL
    ) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        is.data.frame(zscore_table),
        is.character(zscore_columns),
        assertthat::is.string(output_file) || is.null(output_file),
        assertthat::is.string(array_id) || is.null(array_id)
    )

    # Make sure the z-score table has the expected columns
    expected_cols <- c(
        "probe_id",
        "probe_type",
        "probe_sequence",
        "probe_set",
        "snv_position",
        "snv_nucleotide",
        zscore_columns
    )
    if (!all(expected_cols %in% colnames(zscore_table))) {
        stop(
            "zscore_table is missing one or more expected columns\n",
            "Expected columns: ",
            paste(expected_cols, collapse = ", "),
            call. = FALSE
        )
    }

    # If array_id is NULL, switch it to NA_character_
    if (is.null(array_id)) {
        array_id <- NA_character_
    }

    # Make a table of motif data
    motif_table <-
        zscore_table %>%

        # Keep only the seed probe row for each non-background probe set
        dplyr::filter(snv_position == 0) %>%

        # Keep only the probe set name and probe sequence columns
        dplyr::select(probe_set, probe_sequence) %>%

        # Get all possible combinations of probe set names and PBM conditions
        tidyr::expand_grid(pbm_conditions = zscore_columns) %>%

        # Make the z-score motif for each probe set/condition combination
        dplyr::mutate(
            zscore_motif = purrr::map2(
                probe_set,
                pbm_conditions,
                make_zscore_motif,
                zscore_table = zscore_table
            )
        )

    # Convert the data frame into a list of CoRecMotifs
    corecmotifs <-
        purrr::pmap(
            list(
                seed_name = motif_table$probe_set,
                pbm_condition = motif_table$pbm_conditions,
                zscore_motif = motif_table$zscore_motif,
                seed_sequence = motif_table$probe_sequence,
                array_id = array_id
            ),
            CoRecMotif
        )

    # Try to save the CoRecMotifs as an RDS file if necessary
    try_catch_save_output(corecmotifs, output_file, "rds")

    # Return the list of CoRecMotifs
    return(corecmotifs)
}

#' Create a z-score motif for a given probe set and PBM condition
#'
#' Creates a z-score motif for a given probe set and PBM condition. The rows
#' correspond to nucleotides and the columns correspond to positions in the
#' motif.
#'
#' @inheritParams zscore_table_to_corecmotifs
#' @param probe_set_name a character string containing the name of the probe set
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
#' print("FILL THIS IN")
make_zscore_motif <- function(zscore_table, probe_set_name, pbm_condition) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        is.data.frame(zscore_table),
        assertthat::is.string(probe_set_name),
        assertthat::is.string(pbm_condition)
    )

    # Make sure the z-score table has the expected columns
    expected_cols <- c(
        "probe_id",
        "probe_type",
        "probe_sequence",
        "probe_set",
        "snv_position",
        "snv_nucleotide",
        pbm_condition
    )
    if (!all(expected_cols %in% colnames(zscore_table))) {
        stop(
            "zscore_table is missing one or more expected columns\n",
            "Expected columns: ",
            paste(expected_cols, collapse = ", "),
            call. = FALSE
        )
    }

    # Get the z-score of the seed probe for this seed_name/pbm_condition combo
    seed_zscore <-
        zscore_table %>%

        # Keep only the seed probe row from the relevant probe set
        dplyr::filter(probe_set == probe_set_name & snv_position == 0) %>%

        # Pull the z-score from the relevant PBM condition column
        dplyr::pull(pbm_condition)

    # Make sure there's exactly one seed probe
    if (length(seed_zscore) != 1) {
        stop(
            "Expected 1 seed probe for the probe set ",
            probe_set_name,
            " but found ",
            length(seed_zscore),
            call. = FALSE
        )
    }

    # Fill in a motif data frame with the seed and SV probe z-scores
    zscore_motif <-
        zscore_table %>%

        # Keep only the SV probe rows from the relevant probe set
        dplyr::filter(probe_set == probe_set_name & snv_position > 0) %>%

        # Keep only the nucleotide, position, and z-score columns
        dplyr::select(
            snv_position,
            snv_nucleotide,
            zscore = !!as.symbol(pbm_condition)
        ) %>%

        # Fill in missing SNV_pos_offset/SNV_nuc combos with the seed z-score
        tidyr::complete(
            snv_position,
            snv_nucleotide,
            fill = list(zscore = seed_zscore)
        ) %>%

        # Reformat so the columns are positions and the rows are nucleotides
        tidyr::pivot_wider(
            names_from = snv_position,
            values_from = zscore
        ) %>%

        # Make the nucleotides the row names
        tibble::column_to_rownames("snv_nucleotide") %>%

        # Convert to a matrix
        as.matrix()

    # Return the z-score motif
    return(zscore_motif)
}

