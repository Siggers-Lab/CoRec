#' Create \linkS4class{CoRecMotif} objects from a table of z-score data
#'
#' Creates a list of \linkS4class{CoRecMotif} objects for all possible
#' combinations of the probe sets present in \code{zscore_table} and the PBM
#' conditions given in \code{zscore_columns}.
#'
#' @param zscore_table a data frame of z-scores and annotations for each probe.
#'   See 'Details' of \code{\link{annotate_fluorescence_table}} for a
#'   description of the expected annotation columns.
#' @param zscore_columns a character vector specifying the names of the columns
#'   of \code{zscore_table} that contain z-score data.
#' @param output_file the path to the RDS file where the list of CoRecMotif
#'   objects will be written. If NULL (the default), no file is written.
#' @param array_id an optional (but recommended) tag specifying the particular
#'   array/experiment the fluorescence data is from.
#'
#' @return A list of \linkS4class{CoRecMotif} objects, one for each possible
#'   combination of the probe sets in \code{zscore_table} and the PBM conditions
#'   listed in \code{zscore_columns}.
#'
#' @export
#'
#' @examples
#' # Load the example z-score table
#' zscore_table <-
#'     read.table(
#'         "example_data/output/example_rep1_v1_a11_run1_zscores.tsv",
#'         header = TRUE,
#'         sep = "\t",
#'         stringsAsFactors = FALSE
#'     )
#'
#' # Make a list of CoRecMotif objects
#' corecmotifs <-
#'     zscore_table_to_corecmotifs(
#'         zscore_table,
#'         zscore_columns = = c(
#'             "UT_SUDHL4_SWISNF_mix",
#'             "UT_SUDHL4_HDAC1_mix",
#'             "UT_SUDHL4_PRMT5",
#'             "UT_SUDHL4_JMJD2A"
#'         ),
#'         array_id = "v1_a11_run1"
#'     )
zscore_table_to_corecmotifs <-
    function(
        zscore_table,
        zscore_columns,
        output_file = NULL,
        array_id = NULL
    ) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(is.data.frame(zscore_table))
    assertthat::assert_that(is.character(zscore_columns))
    assertthat::assert_that(
        assertthat::is.string(output_file) || is.null(output_file),
        msg = "output_file is not a character vector or NULL"
    )
    assertthat::assert_that(
        assertthat::is.string(array_id) || is.null(array_id),
        msg = "array_id is not a character vector or NULL"
    )

    # Make sure the z-score table has the expected columns
    expected_cols <- c(
        "probeID",
        "probe_type",
        "probe_seq",
        "seed_names",
        "SNV_pos_offset",
        "SNV_nuc",
        zscore_columns
    )
    if (! all(expected_cols %in% colnames(zscore_table))) {
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

    # Convert the data frame into a list of CoRecMotif objects
    corecmotifs <-
        purrr::pmap(
            list(
                seed_name = motif_table$seed_names,
                pbm_condition = motif_table$pbm_conditions,
                zscore_motif = motif_table$zscore_motif,
                seed_sequence = motif_table$probe_seq,
                array_id = array_id
            ),
            CoRecMotif
        )

    # Save the list of all CoRecMotifs as an RDS file if necessary
    if (!is.null(output_file)) {
        tryCatch(
            # Try to save the CoRecMotifs to the output file
            suppressWarnings(
                saveRDS(
                    corecmotifs,
                    output_file
                )
            ),
            # If it fails, skip the output saving step with a warning
            error = function(e) {
                warning(
                    "Could not write to output file '",
                    output_file,
                    "'\nSkipping output file creation...",
                    call. = FALSE
                )
            }
        )
    }

    # Return the list of CoRecMotifs
    return(corecmotifs)
}

#' Create a z-score motif for a given probe set and PBM condition
#'
#' Creates a z-score motif for a given probe set and PBM condition. The rows
#' correspond to nucleotides and the columns correspond to positions in the
#' motif.
#'
#' @param zscore_table a data frame of z-scores and annotations for each probe.
#'   See 'Details' of \code{\link{annotate_fluorescence_table}} for a
#'   description of the expected annotation columns.
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
#'         "example_data/output/example_rep1_v1_a11_run1_zscores.tsv",
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
#'         pbm_condition = "UT_SUDHL4_PRMT5"
#'     )
make_zscore_motif <- function(zscore_table, probe_set, pbm_condition) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(is.data.frame(zscore_table))
    assertthat::assert_that(assertthat::is.string(probe_set))
    assertthat::assert_that(assertthat::is.string(pbm_condition))

    # Make sure the z-score table has the expected columns
    expected_cols <- c(
        "probeID",
        "probe_type",
        "probe_seq",
        "seed_names",
        "SNV_pos_offset",
        "SNV_nuc",
        pbm_condition
    )
    if (! all(expected_cols %in% colnames(zscore_table))) {
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
        dplyr::filter(seed_names == probe_set & SNV_pos_offset == 0) %>%

        # Pull the z-score from the relevant PBM condition column
        dplyr::pull(pbm_condition)

    # Make sure there's exactly one seed probe
    if (length(seed_zscore) != 1) {
        stop(
            "Expected 1 seed probe for the probe set ",
            probe_set,
            " but found ",
            length(seed_zscore),
            call. = FALSE
        )
    }

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
        tibble::column_to_rownames("SNV_nuc") %>%

        # Convert to a matrix
        as.matrix()

    # Return the z-score motif
    return(zscore_motif)
}

