#' Transform fluorescence matrix to z-scores
#'
#' Converts a data frame of fluorescence values into PBM condition-wise z-scores
#' using the background fluorescence values.
#'
#' This function normalizes raw fluorescence values against the fluorescence
#' values of the background probes. It first log-transforms all the fluorescence
#' values for a given PBM condition, then calculates a z-score like statistic
#' using the mean and standard deviation of the background probes.
#'
#' @param annotated_fluorescence_matrix A data frame of fluorescence values and
#'   annotations for each probe, such as output by
#'   annotate_fluorescence_matrix().
#'
#' @return A data frame of PBM condition-wise fluorescence value z-scores.
#' @export
#'
#' @examples
#' zscore_matrix <- fluorescence_to_zscore_matrix(annotated_fluorescence_matrix)
fluorescence_to_zscore_matrix <-
    function(
        annotated_fluorescence_matrix
    ) {
        # Find the indices of the rows that contain background probes
        background_rows <-
            which(annotated_fluorescence_matrix$probe_type == "BACKGROUND")

        # Calculate the column-wise z-scores for each PBM condition
        zscore_matrix <-
            annotated_fluorescence_matrix %>%

            # Change only the data columns and not the annotation columns
            dplyr::mutate(dplyr::across(
                dplyr::matches("o1_|o2_|or_|br_"),
                fluorescence_to_zscore_helper,
                background_rows
            ))

        # Return the matrix of z-scores
        return(zscore_matrix)
    }


# Helper function for fluorescence_to_zscore_matrix()
#
# Converts a vector of fluorescence values into z-scores using the background
# fluorescence values.
#
# @param fluorescence A vector of fluorescence values. Generally this will be an
#   individual column representing a single PBM condition from a data frame of
#   fluorescence values.
# @param background_rows A vector of the indices of background probe
#   fluorescence values in \code{fluorescence}.
#
# @return A vector of fluorescence value z-scores.
fluorescence_to_zscore_helper <-
    function(
        fluorescence,
        background_rows
    ) {
        # Replace negative fluorescence values with NA
        fluorescence <-
            replace(fluorescence, which(fluorescence < 0), NA)

        # Log transform the fluorescence values
        log_fluorescence <-
            log(fluorescence)

        # Get the mean of the background probe fluorescence values
        background_mean <-
            mean(log_fluorescence[background_rows], na.rm = TRUE)

        # Get the standard deviation of the background probe fluorescence values
        background_sd <-
            sd(log_fluorescence[background_rows], na.rm = TRUE)

        # Calculate the z-scores of the vector of log(fluorescence values)
        zscores <-
            (log_fluorescence - background_mean) / background_sd

        # Return the vector of z-scores
        return(zscores)
    }

