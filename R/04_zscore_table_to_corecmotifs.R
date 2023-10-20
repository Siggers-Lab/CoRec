#' Create CoRecMotifs from a table of z-score data
#'
#' Converts a table of fluorescence z-scores into a list of
#' [CoRecMotifs][CoRecMotif-class].
#'
#' This function creates a [CoRecMotif][CoRecMotif-class] for every possible
#' combination of the z-score columns listed in `zscore_columns` and the
#' non-background probe sets that appear in the `probe_set` column of
#' `zscore_table`. Each column listed in `zscore_columns` should contain data
#' from one PBM condition, and the column name will be used to populate the
#' `pbm_condition` slot of the resulting [CoRecMotifs][CoRecMotif-class]. The
#' names of the non-background probe sets found in the `probe_set` column will
#' populate the `probe_set` slot, and `array_id` (if provided) will populate the
#' `array_id` slot. If no array ID is provided, a random one will be generated
#' in the format "random_id_xxxxxxxx" where "xxxxxxxx" is replaced with a
#' randomly generated 8 digit number. The `probe_set`, `pbm_condition`, and
#' `array_id` slots will be pasted together to create a uniquely identifying
#' name for each [CoRecMotif][CoRecMotif-class].
#'
#' @inheritParams annotate_fluorescence_table
#' @inheritParams CoRecMotif
#' @param zscore_table `data.frame`. An annotated table of fluorescence
#'   z-scores. See [hTF_v1_annotation] for expected annotation columns.
#' @param zscore_columns `character`. The names of the columns of `zscore_table`
#'   that contain fluorescence z-scores.
#'
#' @return A list of [CoRecMotifs][CoRecMotif-class], one for each possible
#'   combination of the probe sets in `zscore_table$probe_set` and the PBM
#'   conditions listed in `zscore_columns`.
#'
#' @seealso [hTF_v1_annotation] for a description of the probe annotation
#'   columns.
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
        assertthat::is.string(array_id) || is.null(array_id),
        assertthat::is.string(output_file) || is.null(output_file)
    )

    # Set the expected column names
    expected_cols <- c(
        "probe_id",
        "probe_type",
        "probe_sequence",
        "probe_set",
        "snv_position",
        "snv_nucleotide",
        zscore_columns
    )

    # Make sure fluorescence_table has the expected columns and remove extras
    zscore_table <- check_colnames(zscore_table, expected_cols)

    # If no array ID is given, generate a random ID
    # CoRecMotifs have to have different names or check_replicates() won't work
    if (is.null(array_id)) {
        array_id <- create_array_id()
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
                probe_set = motif_table$probe_set,
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

#' Create a z-score motif
#'
#' Creates a z-score motif for a given probe set and PBM condition. The rows
#' correspond to nucleotides and the columns correspond to positions in the
#' motif.
#'
#' @inheritParams zscore_table_to_corecmotifs
#' @param probe_set_name `character(1)`. The name of the probe set to use.
#' @param pbm_condition `character(1)`. The name of the PBM condition to use.
#'
#' @return A matrix with rows corresponding to nucleotides and columns
#'   corresponding to positions in the motif. Each cell is filled with the PBM
#'   condition-specific z-score of the probe with that nucleotide at that
#'   position.
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

    # Set the expected column names
    expected_cols <- c(
        "probe_id",
        "probe_type",
        "probe_sequence",
        "probe_set",
        "snv_position",
        "snv_nucleotide",
        pbm_condition
    )

    # Make sure fluorescence_table has the expected columns and remove extras
    zscore_table <- check_colnames(zscore_table, expected_cols)

    # Get the z-score of the seed probe for this probe_set/pbm_condition combo
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

        # Missing snv_position/snv_nucleotide should have the seed z-score
        tidyr::complete(
            snv_position,
            snv_nucleotide,
            fill = list(zscore = seed_zscore)
        ) %>%

        # Make sure the snv_positions are in the correct order
        # If the snv_position column wasn't numeric, it might be wrong
        dplyr::arrange(
            as.numeric(snv_position), as.character(snv_nucleotide)
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

