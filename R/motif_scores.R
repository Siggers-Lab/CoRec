#' Score a corecmotif
#'
#' Assign a score to a corecmotif object using one of several available scoring
#' metrics. See 'Details' for descriptions of the available metrics.
#'
#' @param corec_motif The corecmotif object to score.
#' @param method The scoring metric to use. Must be one of ..............
#' @param window_length The length of the sliding window to use for metrics
#'   involving rolling means or medians.
#'
#' @return A corecmotif object with the motif_score and motif_score_type slots
#'   filled in. If the corecmotif object supplied to the \code{corec_motif}
#'   argument already has values in these slots, they will be overwritten.
#' @export
#'
#' @examples
#' scored_corec_motif <- score_motif(corec_motif)
score_motif <-
    function(
        corec_motif,
        method = c(
            "seed_zscore",
            "mean_highest_delta_zscore",
            "median_highest_delta_zscore",
            "rolling_mean_highest_delta_zscore",
            "rolling_median_highest_delta_zscore",
            "mean_information_content",
            "median_information_content",
            "rolling_mean_information_content",
            "rolling_median_information_content"
        ),
        window_length = 5
    ) {
        # Make sure the selected scoring method is a valid option
        method <-
            match.arg(method)

        # Calculate the appropriate motif score
        score <-
            switch(
                method,
                seed_zscore =
                    seed_zscore(corec_motif@zscore_motif),
                mean_highest_delta_zscore = NA,
                median_highest_delta_zscore = NA,
                rolling_mean_highest_delta_zscore = NA,
                rolling_median_highest_delta_zscore = NA,
                mean_information_content = NA,
                median_information_content = NA,
                rolling_mean_information_content = NA,
                rolling_median_information_content = NA
            )

        # Modify the corecmotif object to include the score and score type
        corec_motif@motif_score <- score
        corec_motif@motif_score_type <- method

        # Return the modified corecmotif object
        return(corec_motif)
    }


# Seed probe z-score
#
# Finds the fluorescence value z-score of the seed probe of the z-score motif.
#
# @param zscore_motif A data frame representing a z-score motif, where the rows
#   are nucleotides and the columns are positions in the motif.
#
# @return The fluorescence value z-score of the seed probe of the z-score motif.
seed_zscore <- function(zscore_motif) {
    # The seed probe z-score shows up at every position of the z-score motif
    seed_probe_zscore <-
        Reduce(dplyr::intersect, as.list(zscore_motif))

    # Return the z-score
    return(seed_probe_zscore)
}


# Title
#
# @param ppm
# @param width
#
# @return
rolling_ic <- function(ppm, width = 5) {
    # Convert the PPM to an information content matrix
    icm <-
        universalmotif::convert_type(ppm, type = "ICM")

    #
    max_sliding_window_ic <-
        icm["motif"] %>%
        colSums() %>%
        zoo::rollmean(width) %>%
        max()

    return(max_sliding_window_ic)
}

