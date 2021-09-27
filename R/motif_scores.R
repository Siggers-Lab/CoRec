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


# Rolling mean highest (delta) z-score
#
# Description
#
# @param zscore_motif A data frame representing a z-score motif, where the rows
#   are nucleotides and the columns are positions in the motif.
#
# @return
rolling_mean_highest_zscore <- function(zscore_motif, width = 5) {
    # Calculate the average of the highest scoring probe at each position
    mean_highest_zscore <-
        zscore_motif %>%

        # Find the highest z-score at each position
        dplyr::summarise_all(max) %>%

        # Convert the data frame of highest z-scores to a numeric vector
        as.numeric() %>%

        # Calculate the average of the highest z-scores across a sliding window
        zoo::rollmean(width) %>%

        # Find the maximum value from the vector of rolling averages
        max()

    # Return the average of the highest scoring probe at each position
    return(mean_highest_zscore)
}


# Title
#
# @param zscore_motif
#
# @return
rolling_median_highest_zscore <- function(zscore_motif, width = 5) {
    # Calculate the average of the highest scoring probe at each position
    median_highest_zscore <-
        zscore_motif %>%

        # Find the highest z-score at each position
        dplyr::summarise_all(max) %>%

        # Convert the data frame of highest z-scores to a numeric vector
        as.numeric() %>%

        # Calculate the average of the highest z-scores across a sliding window
        zoo::rollmedian(width) %>%

        # Find the maximum value from the vector of rolling averages
        max()

    # Return the average of the highest scoring probe at each position
    return(median_highest_zscore)
}


# Title
#
# @param zscore_motif
#
# @return
median_highest_zscore <- function(zscore_motif) {
    # Calculate the average of the highest scoring probe at each position
    median_highest_zscore <-
        zscore_motif %>%

        # Find the highest z-score at each position
        dplyr::summarise_all(max) %>%

        # Convert the data frame of highest z-scores to a numeric vector
        as.numeric() %>%

        # Calculate the average of the highest z-scores across a sliding window
        median()

    # Return the average of the highest scoring probe at each position
    return(median_highest_zscore)
}


count_above_threshold <- function(zscore_motif, threshold) {
    count_above_threshold <-
        zscore_motif %>%

        # Find the highest z-score at each position
        dplyr::summarise_all(max) %>%

        # Convert the data frame of highest z-scores to a numeric vector
        as.numeric() %>%

        # Keep only those z-scores that are greater than the threshold
        subset(. > threshold) %>%

        # Count the number of remaining z-scores
        length()

    return(count_above_threshold)
}


# Title
#
# @param zscore_motif
# @param threshold
#
# @return
longest_window <- function(zscore_motif, threshold) {
    # Find the lengths of windows with z-scores greater than the threshold
    window_lengths <-
        zscore_motif %>%

        # Find the highest z-score at each position
        dplyr::summarise_all(max) %>%

        # Convert to logicals with TRUE for z-scores above the threshold
        dplyr::mutate(across(.fns = ~. > threshold)) %>% as.logical() %>%

        # Find the length of each sequence of TRUE and FALSE
        rle()

    # Find the longest sequence of TRUE, or 0 if no z-scores were high enough
    longest_window <-
        ifelse(
            any(window_lengths$values),
            max(window_lengths$lengths[window_lengths$values]),
            0
        )

    return(longest_window)
}

# Title
#
# @param ppm
#
# @return
mean_information_content <- function(ppm) {
    # Average information content = total information content / motif length
    mean_information_content <-
        ppm["icscore"] / ncol(ppm["motif"])

    # Return the mean information content
    return(mean_information_content)
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

# Title
#
# @param ppm
# @param threshold
#
# @return
count_above_ic_threshold <- function(ppm, threshold) {
    icm <-
        universalmotif::convert_type(ppm, type = "ICM")

    count_above_threshold <-
        icm["motif"] %>%

        # Find the total information content at each position
        colSums() %>%

        # Keep only those scores that are greater than the threshold
        subset(. > threshold) %>%

        # Count the number of remaining scores
        length()

    return(count_above_threshold)
}


# Title
#
# @param ppm
# @param threshold
#
# @return
longest_ic_window <- function(ppm, threshold) {
    icm <-
        universalmotif::convert_type(ppm, type = "ICM")

    # Find the lengths of windows with information content above the threshold
    window_lengths <-
        icm["motif"] %>%

        # Find the total information content at each position
        colSums() %>% as.data.frame() %>% t() %>% as.data.frame() %>%

        # Convert to logicals with TRUE for z-scores above the threshold
        dplyr::mutate(across(.fns = ~. > threshold)) %>% as.logical() %>%

        # Find the length of each sequence of TRUE and FALSE
        rle()

    # Find the longest sequence of TRUE, or 0 if no z-scores were high enough
    longest_window <-
        ifelse(
            any(window_lengths$values),
            max(window_lengths$lengths[window_lengths$values]),
            0
        )

    return(longest_window)
}

