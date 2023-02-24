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
#' @inheritParams annotate_fluorescence_table
#' @param fluorescence_table `data.frame`. An annotated table of fluorescence
#'   values. See [hTF_v1_annotation] for expected annotation columns.
#'
#' @return A data frame of fluorescence z-scores and the corresponding probe
#'   information. See [hTF_v1_annotation] for a description of the probe
#'   annotation columns.
#'
#' @seealso [hTF_v1_annotation] for a description of the probe annotation
#'   columns.
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

    # Set the expected column names
    expected_cols <- c(
        "probe_id",
        "probe_type",
        "probe_sequence",
        "probe_set",
        "snv_position",
        "snv_nucleotide",
        fluorescence_columns
    )

    # Make sure fluorescence_table has the expected columns and remove extras
    fluorescence_table <- check_colnames(fluorescence_table, expected_cols)

    # Find the indices of the rows that contain background probes
    background_rows <- which(fluorescence_table$probe_type == "BACKGROUND")

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

