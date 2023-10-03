#' Import fluorescence data
#'
#' Reads a file containing fluorescence data (in which rows correspond to probes
#' and columns correspond to fluorescence values) and reformats it. The output
#' from this function is suitable for use by [annotate_fluorescence_table].
#'
#' The fluorescence data file supplied to `fluorescence_file` must contain
#' \eqn{n + 3} tab-delimited columns, where n is the number of PBM conditions
#' supplied to `pbm_conditions`. This file must NOT have a header. The columns
#' must be as follows:
#'
#' 1. The probe ID. The orientation tags (i.e., "_o1" or "_o2") will be removed
#' from these IDs.
#' 2. The probe sequence. This column will be removed.
#' 3. The number of PBM conditions profiled on this array. This column will be
#' removed.
#'
#' Columns 4 through \eqn{(n + 3)} must contain the fluorescence values for the
#' PBM conditions in the order they appear in `pbm_conditions`.
#'
#' @param fluorescence_file `character(1)`. The path to the fluorescence data
#'   file to load. See 'Details' for expected columns.
#' @param pbm_conditions `character`. The PBM conditions (e.g., cell type,
#'   treatment, and factor profiled) in the order they appear in
#'   `fluorescence_file`.
#'
#' @return A data frame containing the fluorescence values from the specified
#'   fluorescence data file. The first column will be named "probe_id", and the
#'   remaining columns will be named according to the PBM conditions supplied to
#'   `pbm_conditions`.
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN")
load_fluorescence_data <- function(fluorescence_file, pbm_conditions) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        assertthat::is.string(fluorescence_file) &&
            file.exists(fluorescence_file),
        is.character(pbm_conditions)
    )

    # Make sure the PBM conditions are valid column names
    pbm_conditions_fixed <- make.names(pbm_conditions)

    # Throw a warning if the column names were changed
    if (!identical(pbm_conditions, pbm_conditions_fixed)) {
        warning(
            "pbm_conditions contained invalid column name(s)\n",
            "New names:",
            lapply(1:length(pbm_conditions), function(i) {
                if (pbm_conditions[i] != pbm_conditions_fixed[i]) {
                    paste0(
                        "\n\t",
                        pbm_conditions[i],
                        " -> ",
                        pbm_conditions_fixed[i]
                    )
                }
            }),
            call. = FALSE
        )
    }

    # Load the table of fluorescence values
    fluorescence_table <-
        utils::read.table(
            fluorescence_file,
            header = FALSE,
            sep = "\t",
            strip.white = TRUE,
            stringsAsFactors = FALSE
        )

    # Make sure the list of PBM conditions is the right length
    if (ncol(fluorescence_table) != (length(pbm_conditions_fixed) + 3)) {
        stop(
            "pbm_conditions is the wrong length\n",
            "Expected ",
            ncol(fluorescence_table) -3,
            " values but got ",
            length(pbm_conditions_fixed),
            call. = FALSE
        )
    }

    # Reformat the fluorescence table
    fluorescence_table <-
        fluorescence_table %>%

        # Rename the columns
        magrittr::set_colnames(
            c(
                "probe_id",
                "probe_sequence",
                "n_conditions",
                pbm_conditions_fixed
            )
        ) %>%

        # Remove the unnecessary columns
        dplyr::select(-probe_sequence, -n_conditions) %>%

        # Remove the orientation tag from the probe IDs
        dplyr::mutate(probe_id = gsub("_o[1-2]", "", probe_id))

    # Return the fluorescence table
    return(fluorescence_table)
}

