#' Transform fluorescence values into z-scores
#'
#' Converts an annotated data frame of fluorescence data into z-scores.
#'
#' Raw fluorescence values are normalized against the fluorescence values of the
#' background probes. Specifically, all the fluorescence values for a given PBM
#' condition are log-transformed. Then a z-score-like statistic is calculated
#' using the mean and standard deviation of the fluorescence values of the
#' background probes in that condition.
#'
#' @param fluorescence_table a data frame of fluorescence values and annotations
#'   for each probe. See 'Details' of [annotate_fluorescence_table()] for a
#'   description of the expected annotation columns.
#' @param fluorescence_columns a character vector specifying the names of the
#'   columns of `fluorescence_table` that contain fluorescence data.
#' @param output_file the path to the TSV file where the annotated z-score table
#'   will be written. If NULL, no file is written. (Default: NULL)
#'
#' @return A data frame of column-wise fluorescence z-scores and probe
#'   annotations.
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN")
fluorescence_to_zscore_table <-
    function(
        fluorescence_table,
        fluorescence_columns,
        output_file = NULL
    ) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        is.data.frame(fluorescence_table),
        is.character(fluorescence_columns),
        assertthat::is.string(output_file) || is.null(output_file)
    )

    # Make sure the fluorescence table has the expected columns
    expected_cols <- c(
        "probe_id",
        "probe_type",
        "probe_sequence",
        "probe_set",
        "snv_position",
        "snv_nucleotide",
        fluorescence_columns
    )
    if (!all(expected_cols %in% colnames(fluorescence_table))) {
        stop(
            "fluorescence_table is missing one or more expected columns\n",
            "Expected columns: ",
            paste(expected_cols, collapse = ", "),
            call. = FALSE
        )
    }

    # Find the indices of the rows that contain background probes
    background_rows <-
        which(fluorescence_table$probe_type == "BACKGROUND")

    # Calculate the column-wise z-scores for each PBM condition
    zscore_table <-
        fluorescence_table %>%

        # Change only the data columns and not the annotation columns
        dplyr::mutate(
            # Replace negative fluorescence values with NA
            dplyr::across(
                dplyr::matches(fluorescence_columns),
                ~ replace(.x, which(.x < 0), NA)
            ),
            # Replace NAs with the minimum fluorescence value
            # This will catch anything that was negative or NA
            dplyr::across(
                dplyr::matches(fluorescence_columns),
                ~ replace(.x, which(is.na(.x)), min(.x, na.rm = TRUE))
            ),
            # Log transform the fluorescence values
            dplyr::across(
                dplyr::matches(fluorescence_columns),
                log
            ),
            # Calculate the z-scores of the log(fluorescence values)
            dplyr::across(
                dplyr::matches(fluorescence_columns),
                ~ (.x - mean(.x[background_rows])) / sd(.x[background_rows])
            )
        )

    # Try to save the z-score table if necessary
    try_catch_save_output(zscore_table, output_file, "tsv")

    # Return the table of z-scores
    return(zscore_table)
}

