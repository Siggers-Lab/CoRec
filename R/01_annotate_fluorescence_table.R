#' Annotate fluorescence data
#'
#' Adds probe information to a table of fluorescence data from a CoRec
#' experiment.
#'
#' The fluorescence data table must have a column named "`probe_id`". The
#' remaining column names must match the names listed in `fluorescence_columns`.
#' Any additional columns will be dropped.
#'
#' See [hTF_v1_annotation] for a description of the expected columns of the
#' annotation table. Any additional columns will be dropped.
#'
#' The IDs in the "`probe_id`" column in the annotation table must match the IDs
#' in the "`probe_id`" column of the fluorescence table. IDs that are present in
#' the annotation table but not the fluorescence table will be dropped silently.
#' IDs that are present in the fluorescence table but not the annotation table
#' will be dropped with a warning.
#'
#' @param fluorescence_table `data.frame`. A table of fluorescence values.
#' @param fluorescence_columns `character`. The names of the columns of
#'   `fluorescence_table` that contain fluorescence data.
#' @param annotation `data.frame`. The probe annotations to add to the
#'   fluorescence table. See [hTF_v1_annotation] for expected annotation
#'   columns.
#' @param output_file `character(1)` or `NULL`. The path to the file where the
#'   results should be saved, or NULL not to save the results. (Default: NULL)
#'
#' @seealso [hTF_v1_annotation] for a description of the probe annotation
#'   columns.
#'
#' @return A data frame of fluorescence values and the corresponding probe
#'   information. See [hTF_v1_annotation] for a description of the probe
#'   annotation columns.
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN")
annotate_fluorescence_table <-
    function(
        fluorescence_table,
        fluorescence_columns,
        annotation,
        output_file = NULL
    ) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        is.data.frame(fluorescence_table),
        is.character(fluorescence_columns),
        is.data.frame(annotation),
        assertthat::is.string(output_file) || is.null(output_file)
    )

    # Set the expected column names
    expected_fluo_cols <- c("probe_id", fluorescence_columns)
    expected_anno_cols <- c(
        "probe_id",
        "probe_type",
        "probe_sequence",
        "probe_set",
        "snv_position",
        "snv_nucleotide"
    )

    # Make sure fluorescence_table has the expected columns and remove extras
    fluorescence_table <- check_colnames(fluorescence_table, expected_fluo_cols)

    # Make sure annotation has the expected columns and remove extras
    annotation <- check_colnames(annotation, expected_anno_cols)

    # Check how many rows the fluorescence table has
    n_rows <- nrow(fluorescence_table)

    # Add the probe annotations to the fluorescence table
    fluorescence_table <-
        dplyr::inner_join(annotation, fluorescence_table, by = c("probe_id"))

    # Give a warning if the fluorescence table lost any rows
    if (nrow(fluorescence_table) < n_rows) {
        warning(
            "annotation is missing probe IDs present in fluorescence_table\n",
            "Are you sure you're using the correct annotation table?",
            call. = FALSE
        )
    }

    # Try to save the annotated fluorescence table if necessary
    try_catch_save_output(fluorescence_table, output_file, "tsv")

    # Return the annotated fluorescence table
    return(fluorescence_table)
}

