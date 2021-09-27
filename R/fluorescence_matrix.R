#' Import fluorescence matrix
#'
#' Reads a set of individual o1, o2, or, and br fluorescence matrices and
#' creates a single data frame from them, with rows corresponding to probes and
#' columns corresponding to PBM conditions.
#'
#' This function takes as input the fluorescence matrices output by the PBM
#' preprocessing pipeline. It expects the filenames to be of the form
#' \code{<matrix_base_name>_<probe_orientation>.dat}, where
#' \code{<matrix_base_name>} is the value passed to the \code{matrix_base_name}
#' parameter and \code{<probe_orientation>} is \verb{"o1"}, \verb{"o2"},
#' \verb{"or"}, or \verb{"br"} (all four files must be present). These files
#' should contain \eqn{3 + n} tab-delimited columns, where n is the number of
#' PBM conditions supplied to the \code{pbm_conditions} parameter. The first
#' column must contain the probe IDs, and the fourth through nth columns must
#' contain the fluorescence values for the PBM conditions in the order listed.
#' The second and third columns are removed. The files should not contain
#' headers.
#'
#' The output data frame contains \eqn{1 + (4 * n)} columns, where n is again
#' the number of PBM conditions supplied to the \code{pbm_conditions} parameter.
#' The first column will be named \verb{"probeID"} and will contain the probe
#' IDs with their orientation tags (i.e., \verb{"_o1"} or \verb{"_o2"}) removed.
#' The remaining columns will contain the fluorescence values for each probe
#' orientation and PBM condition and will be named based on the vector of PBM
#' conditions supplied to \code{pbm_conditions} and the tag supplied to
#' \code{run_tag} (if any). If a tag is provided, the column names will be of
#' the form \verb{"<run_tag>_<probe_orientation>_<pbm_condition>"}. If no tag is
#' provided, the column names will be of the form
#' \verb{"<probe_orientation>_<pbm_condition>"}.
#'
#' @param matrix_directory The path to the directory containing the individual
#'   o1, o2, or, and br fluorescence matrices to import.
#' @param matrix_base_name The base filename of the fluorescence matrices to
#'   import.
#' @param pbm_conditions A character vector of the names of the PBM conditions.
#'   These will be used as column names.
#' @param run_tag An optional (but recommended) tag specifying the particular
#'   experiment these fluorescence matrices are from. It will be incorporated
#'   into the column names along with the PBM conditions. This option is useful
#'   for differentiating replicates when merging two matrices downstream.
#'
#' @return A data frame containing the fluorescence values from the specified
#'   set of o1, o2, or, and br matrices.
#' @export
#'
#' @examples
#' fluorescence_matrix <- read_fluorescence_matrices(
#'     matrix_directory = paste0(
#'         "/projectnb/siggers/data/hTF_array_project/data/data_matrices/",
#'         "v1_a11_run1"
#'     ),
#'     matrix_base_name = "hTF_v1_SUDHL4_14jan21",
#'     pbm_conditions = c(
#'         "UT_SUDHL4_SMARCA4MIX",
#'         "UT_SUDHL4_HDAC1MIX",
#'         "UT_SUDHL4_SUZ12",
#'         "UT_SUDHL4_PRMT5",
#'         "UT_SUDHL4_JMJD2A",
#'         "UT_SUDHL4_BMI1",
#'         "UT_SUDHL4_DNMT3A",
#'         "UT_SUDHL4_ASH2L"
#'     )
#' )
read_fluorescence_matrices <-
    function(
        matrix_directory,
        matrix_base_name,
        pbm_conditions,
        run_tag = NA
    ) {
        # Load the separate matrices for each probe orientation into a list
        fluorescence_matrices_list <-
            lapply(c("o1", "o2", "or", "br"), function(probe_orientation) {
                # Put together the file name of the matrix for this orientation
                file_name <-
                    paste0(
                        matrix_directory,
                        "/",
                        matrix_base_name,
                        "_",
                        probe_orientation,
                        ".dat"
                    )

                # Make column names based on the orientation and PBM conditions
                pbm_conditions_columns <-
                    paste(probe_orientation, pbm_conditions, sep = "_")

                # If a run tag is provided, add it to the column names
                if (!is.na(run_tag)) {
                    pbm_conditions_columns <-
                        paste(run_tag, pbm_conditions_columns, sep = "_")
                }

                # Load the matrix
                data_matrix <-
                    read.table(
                        file_name,
                        header = FALSE,
                        sep = "\t",
                        strip.white = TRUE,
                        stringsAsFactors = FALSE
                    ) %>%

                    # Rename the columns
                    magrittr::set_colnames(
                        c(
                            "probeID",
                            "dummy1",
                            "dummy2",
                            pbm_conditions_columns
                        )
                    ) %>%

                    # Remove the unnecessary columns
                    dplyr::select(-dummy1,-dummy2) %>%

                    # Remove the orientation tag from the probe IDs
                    dplyr::mutate(probeID = gsub("_o[1-2]", "", probeID))

            # Return the matrix for this probe orientation
            return(data_matrix)
        })

        # Merge the list of fluorescence matrices into a single matrix
        fluorescence_matrix <-
            plyr::join_all(fluorescence_matrices_list, by = "probeID")

        # Return the full combined fluorescence matrix
        return(fluorescence_matrix)
    }


#' Annotate fluorescence matrix
#'
#' Adds probe information to a data frame of fluorescence values. See 'Details'
#' for a description of the inputs, including required annotation columns.
#'
#' The input fluorescence matrix is expected to be of the form output by
#' read_fluorescence_matrices(), with rows corresponding to probes and columns
#' corresponding to PBM conditions. See \code{\link{read_fluorescence_matrices}}
#' for a more detailed description of the expected format.
#'
#' The input annotation table must contain the columns \code{probeID},
#' \code{probe_or}, \code{probe_repl}, \code{probe_type}, ... The expected
#' content of these columns is described in detail below.
#'
#' \describe{
#'
#' \item{\code{probeID}}{The probe ID. This column will be used to map the probe
#' annotations to the values in the fluorescence matrix, so the IDs in this
#' column must match the IDs in the probeID column in the fluorescence matrix.
#' IDs that are not present in both will be dropped.}
#'
#' \item{\code{probe_or}}{The orientation of the probe relative to the slide.
#' This column should contain the values "o1" and "o2". Only "o1" rows will be
#' used.}
#'
#' \item{\code{probe_repl}}{The replicate number of the probe. This column
#' should contain the values "r1", "r2", etc. Only "r1" rows will be used.}
#'
#' \item{\code{probe_type}}{The probe type. This column should contain the
#' values "MOTIF" and "BACKGROUND".}
#'
#' }
#'
#' @param annotation A data frame of PBM probe annotations. See 'Details' for
#'   expected columns.
#' @param fluorescence_matrix A data frame of fluorescence values, such as
#'   output by read_fluorescence_matrices().
#'
#' @return A data frame of fluorescence values and annotations for each probe.
#' @export
#'
#' @examples
#' annotated_fluorescence_matrix <- annotate_fluorescence_matrix(
#'     fluorescence_matrix,
#'     annotation
#' )
annotate_fluorescence_matrix <-
    function(
        fluorescence_matrix,
        annotation
    ) {
        # Keep only one annotation row per unique target sequence
        annotation <-
            annotation %>%
            dplyr::filter(probe_repl == "r1" & probe_or == "o1")

        # Merge the annotation table with the fluorescence matrix
        fluorescence_matrix <-
            fluorescence_matrix %>%
            dplyr::inner_join(annotation, ., by = "probeID") %>%

            # Move the probeID column to the beginning of the data frame
            dplyr::relocate(probeID)

        # Return the annotated fluorescence matrix
        return(fluorescence_matrix)
    }

