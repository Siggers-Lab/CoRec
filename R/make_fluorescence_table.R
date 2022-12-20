#' Import and annotate fluorescence data
#'
#' Reads a file containing fluorescence data (in which rows correspond to probes
#' and columns correspond to fluorescence values) and adds probe metadata from
#' an annotation file.
#'
#' ## Fluorescence File
#' The fluorescence data file supplied to the \code{fluorescence_file} parameter
#' must contain \eqn{3 + n} tab-delimited columns, where n is the number of PBM
#' conditions supplied to the \code{pbm_conditions} parameter. This file must
#' NOT have a header. The columns must be as follows:
#'
#' \describe{
#'
#' \item{\code{Column 1}}{The probe ID. The orientation tags (i.e., \verb{"_o1"}
#' or \verb{"_o2"}) will be removed from these IDs.}
#'
#' \item{\code{Column 2}}{The probe sequence.}
#'
#' \item{\code{Column 3}}{The number of PBM conditions profiled on this array.}
#'
#' \item{\code{Columns 4 through \eqn{(n+3)}}}{The fluorescence values for the
#' PBM conditions in the order they appear in the \code{pbm_conditions}
#' parameter.}
#'
#' }
#'
#' ## Annotation File
#' The probe annotation file supplied to the \code{annotation_file} parameter
#' must have a header and must contain tab-delimited columns named
#' \code{probeID}, \code{probe_type}, \code{probe_seq}, \code{seed_names},
#' \code{SNV_pos_offset}, and \code{SNV_nuc}. The expected content of these
#' columns is described below. The ordering of the columns does not matter. If
#' additional columns are present, they will be removed.
#'
#' \describe{
#'
#' \item{\code{probeID}}{The probe ID.}
#'
#' \item{\code{probe_type}}{The probe type. This column should contain the
#' values "MOTIF" and "BACKGROUND".}
#'
#' \item{\code{probe_seq}}{The probe sequence.}
#'
#' \item{\code{seed_names}}{The name of the probe set this probe is part of.}
#'
#' \item{\code{SNV_pos_offset}}{The position of the SNV in the sequence.}
#'
#' \item{\code{SNV_nuc}}{The nucleotide at the SNV position.}
#'
#' }
#'
#' The IDs in the \code{probeID} column must match the IDs in the first column
#' of the fluorescence data file (after the orientation tag has been removed).
#' IDs that are not present in both will be dropped.The output data frame
#' contains \eqn{6 + n} columns, where n is again the number of PBM conditions
#' supplied to the \code{pbm_conditions} parameter. The first six columns will
#' be the columns from the annotation file (described above). The remaining
#' columns will contain the fluorescence values for each condition and will be
#' named based on the vector of PBM conditions supplied to \code{pbm_conditions}
#' and the tag supplied to \code{array_id} (if any).
#'
#' @param fluorescence_file the name of the file containing the fluorescence
#'   data to load. See 'Details' for expected columns.
#' @param pbm_conditions a character vector specifying the PBM conditions (e.g.,
#'   cell type, treatment, and factor profiled) in the order they appear in
#'   \code{fluorescence_file}.
#' @param annotation_file the name of the file containing the probe annotations
#'   to use. See 'Details' for expected columns.
#' @param array_id an optional (but recommended) tag specifying the particular
#'   array/experiment the fluorescence data is from.
#'
#' @return A data frame containing the fluorescence values from the specified
#'   fluorescence data file and the probe annotations from the specified
#'   annotation file. See 'Details' for a description of each column.
#'
#' @export
#'
#' @examples
#' # Load and annotate the example fluorescence table
#' fluorescence_table <-
#'     make_fluorescence_table(
#'         fluorescence_file = "example_data/hTF_v1_example_fluorescence.dat",
#'         pbm_conditions = c(
#'             "UT_SUDHL4_SMARCA4MIX",
#'             "UT_SUDHL4_HDAC1MIX",
#'             "UT_SUDHL4_SUZ12",
#'             "UT_SUDHL4_PRMT5"
#'         ),
#'         annotation_file = "example_data/hTF_v1_annotation.tsv",
#'         array_id = "v1_a6_run1"
#'     )
make_fluorescence_table <-
    function(
        fluorescence_file,
        pbm_conditions,
        annotation_file,
        array_id = NA
    ) {
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

    # Add the run tag (if provided) to the PBM conditions to make column names
    if (is.na(array_id)) {
        column_names <- pbm_conditions
    } else {
        column_names <- paste(array_id, pbm_conditions, sep = "_")
    }

    # Load the annotation file
    annotation <-
        read.table(
            annotation_file,
            header = TRUE,
            sep = "\t",
            strip.white = TRUE,
            stringsAsFactors = FALSE
        ) %>%

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

    # Reformat the fluorescence table and add the probe annotations
    fluorescence_table <-
        fluorescence_table %>%

        # Rename the columns
        magrittr::set_colnames(
            c(
                "probeID",
                "probe_seq",
                "num_conditions",
                column_names
            )
        ) %>%

        # Remove the unnecessary columns
        dplyr::select(-probe_seq, -num_conditions) %>%

        # Remove the orientation tag from the probe IDs
        dplyr::mutate(probeID = gsub("_o[1-2]", "", probeID)) %>%

        # Merge with the annotation table
        dplyr::inner_join(annotation, ., by = c("probeID"))

    # Return the annotated fluorescence table
    return(fluorescence_table)
}

