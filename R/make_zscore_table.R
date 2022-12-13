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
#'   for each probe. See the return value of
#'   \code{\link{make_fluorescence_table}} for a description of the expected
#'   columns.
#' @param fluorescence_columns a character vector specifying the names of the
#'   columns of \code{fluorescence_table} that contain fluorescence data.
#'
#' @return A data frame of column-wise fluorescence z-scores and probe
#'   annotations.
#' @export
#'
#' @examples
#' fluorescence_table <- make_fluorescence_table(
#'     fluorescence_file =
#'         "example_data/data_matrices/hTF_v1_SUDHL4_14jan21_br.dat",
#'     pbm_conditions = c(
#'         "UT_SUDHL4_SMARCA4MIX",
#'         "UT_SUDHL4_HDAC1MIX",
#'         "UT_SUDHL4_SUZ12",
#'         "UT_SUDHL4_PRMT5",
#'         "UT_SUDHL4_JMJD2A",
#'         "UT_SUDHL4_BMI1",
#'         "UT_SUDHL4_DNMT3A",
#'         "UT_SUDHL4_ASH2L"
#'     ),
#'     annotation_file = "example_data/hTF_v01_PBM_ANNOT.txt",
#'     array_id = "v1_a6_run1"
#' )
#'
#' zscore_table <- make_zscore_table(
#'     fluorescence_table = fluorescence_table,
#'     fluorescence_columns = c(
#'         "v1_a11_run1_UT_SUDHL4_SMARCA4MIX",
#'         "v1_a11_run1_UT_SUDHL4_HDAC1MIX",
#'         "v1_a11_run1_UT_SUDHL4_SUZ12",
#'         "v1_a11_run1_UT_SUDHL4_PRMT5",
#'         "v1_a11_run1_UT_SUDHL4_JMJD2A",
#'         "v1_a11_run1_UT_SUDHL4_BMI1",
#'         "v1_a11_run1_UT_SUDHL4_DNMT3A",
#'         "v1_a11_run1_UT_SUDHL4_ASH2L"
#'     )
#' )
make_zscore_table <-
    function(
        fluorescence_table,
        fluorescence_columns
    ) {
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

        # Return the table of z-scores
        return(zscore_table)
    }

