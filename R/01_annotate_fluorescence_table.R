#' Import and annotate fluorescence data
#'
#' Reads a file containing fluorescence data (in which rows correspond to probes
#' and columns correspond to fluorescence values) and adds probe metadata from
#' an annotation file.
#'
#' ## Fluorescence File
#'
#' The fluorescence data file supplied to `fluorescence_file` must contain
#' \eqn{n + 3} tab-delimited columns, where n is the number of PBM conditions
#' supplied to `pbm_conditions`. This file must NOT have a header. The columns
#' must be as follows:
#'
#' 1. The probe ID. The orientation tags (i.e., "_o1" or "_o2") will be removed
#' from these IDs.
#' 2. The probe sequence.
#' 3. The number of PBM conditions profiled on this array.
#'
#' Columns 4 through \eqn{(n + 3)} must contain the fluorescence values for the
#' PBM conditions in the order they appear in `pbm_conditions`.
#'
#' ## Annotation File
#'
#' The probe annotation file supplied to the `annotation_file` parameter must
#' have a header and must contain tab-delimited columns named `probeID`,
#' `probe_type`, `probe_seq`, `seed_names`, `SNV_pos_offset`, and `SNV_nuc`. The
#' expected content of these columns is described below. The ordering of the
#' columns does not matter. If additional columns are present, they will be
#' removed.
#'
#' * `probeID`: The probe ID.
#' * `probe_type`: The probe type. This column should contain the values "MOTIF"
#'    and "BACKGROUND".
#' * `probe_seq`: The probe sequence
#' * `seed_names`: The seed name (the name of the probe set this probe is part
#'    of).
#' * `SNV_pos_offset`: The position of the SNV in the sequence.
#' * `SNV_nuc`: The nucleotide at the SNV position.
#'
#' The IDs in the `probeID` column must match the IDs in the first column of the
#' fluorescence data file (after the orientation tag has been removed). IDs that
#' are not present in both will be dropped.The output data frame contains \eqn{6
#' + n} columns, where n is again the number of PBM conditions supplied to the
#' `pbm_conditions` parameter. The first six columns will be the columns from
#' the annotation file (described above). The remaining columns will contain the
#' fluorescence values for each condition and will be named based on the vector
#' of PBM conditions supplied to `pbm_conditions`.
#'
#' @param fluorescence_file the path to the file containing the fluorescence
#'   data to load. See 'Details' for expected columns.
#' @param pbm_conditions a character vector specifying the PBM conditions (e.g.,
#'   cell type, treatment, and factor profiled) in the order they appear in
#'   `fluorescence_file`.
#' @param annotation_file the path to the file containing the probe annotations
#'   to use. See 'Details' for expected columns.
#' @param output_file the path to the TSV file where the annotated fluorescence
#'   table will be written. If NULL, no file is created. (Default: NULL)
#'
#' @return A data frame containing the fluorescence values from the specified
#'   fluorescence data file and the probe annotations from the specified
#'   annotation file. See 'Details' for a description of each column.
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN")
annotate_fluorescence_table <-
    function(
        fluorescence_file,
        pbm_conditions,
        annotation_file,
        output_file = NULL
    ) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        assertthat::is.string(fluorescence_file) &&
            file.exists(fluorescence_file),
        is.character(pbm_conditions),
        assertthat::is.string(annotation_file) &&
            file.exists(annotation_file),
        assertthat::is.string(output_file) || is.null(output_file)
    )

    # Load the table of fluorescence values
    fluorescence_table <-
        read.table(
            fluorescence_file,
            header = FALSE,
            sep = "\t",
            strip.white = TRUE,
            stringsAsFactors = FALSE
        )

    # Make sure the list of PBM conditions is the right length
    if (ncol(fluorescence_table) != (length(pbm_conditions) + 3)) {
        stop(
            "pbm_conditions is the wrong length\n",
            "Expected ",
            ncol(fluorescence_table) -3,
            " values but got ",
            length(pbm_conditions),
            call. = FALSE
        )
    }

    # Load the annotation file
    annotation <-
        read.table(
            annotation_file,
            header = TRUE,
            sep = "\t",
            strip.white = TRUE,
            stringsAsFactors = FALSE
        )

    # Make sure the annotation table has the expected column names
    expected_cols <- c(
        "probeID",
        "probe_type",
        "probe_seq",
        "seed_names",
        "SNV_pos_offset",
        "SNV_nuc"
    )
    if (!all(expected_cols %in% colnames(annotation))) {
        stop(
            "annotation_file is missing one or more expected columns\n",
            "Expected columns: ",
            paste(expected_cols, collapse = ", "),
            call. = FALSE
        )
    }

    # Filter the annotation table
    annotation <-
        annotation %>%

        # Select only the required columns
        dplyr::select(
            probeID,
            probe_type,
            probe_seq,
            seed_names,
            SNV_pos_offset,
            SNV_nuc
        ) %>%

        # Group by the probe ID
        dplyr::group_by(probeID) %>%

        # Keep only one row for each probe ID
        dplyr::slice_head() %>%

        # Remove the grouping just in case
        dplyr::ungroup()

    # Check how many rows the fluorescence table has
    n_rows <- nrow(fluorescence_table)

    # Reformat the fluorescence table and add the probe annotations
    fluorescence_table <-
        fluorescence_table %>%

        # Rename the columns
        magrittr::set_colnames(
            c("probeID", "probe_seq", "num_conditions", pbm_conditions)
        ) %>%

        # Remove the unnecessary columns
        dplyr::select(-probe_seq, -num_conditions) %>%

        # Remove the orientation tag from the probe IDs
        dplyr::mutate(probeID = gsub("_o[1-2]", "", probeID)) %>%

        # Merge with the annotation table
        dplyr::inner_join(annotation, ., by = c("probeID"))

    # Give a warning if the fluorescence table lost any rows
    if (nrow(fluorescence_table) < n_rows) {
        warning(
            "Annotation file is missing probeIDs present in fluorescence file",
            "\nAre you sure you're using the correct annotation file?",
            call. = FALSE
        )
    }

    # Try to save the annotated fluorescence table if necessary
    try_catch_save_output(fluorescence_table, output_file, "tsv")

    # Return the annotated fluorescence table
    return(fluorescence_table)
}

