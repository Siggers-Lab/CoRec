shuffle_zscore_motif <- function(zscore_motif) {
    # Figure out the seed probe z-score
    seed_probe_zscore <-
        seed_zscore(zscore_motif)

    # Figure out the seed nucleotide at the first position in the z-score motif
    first_seed_nuc <-
        zscore_motif %>%

        # Keep only the row where the first position has the seed probe z-score
        dplyr::filter(`1` == seed_probe_zscore) %>%

        # Get the name of this row
        rownames()

    shuffled_zscore_motif <-
        zscore_motif %>%

        # Convert the nucleotides into a column instead of rownames
        tibble::rownames_to_column("SNV_nuc") %>%

        # Convert to long format
        tidyr::pivot_longer(where(is.numeric)) %>%

        # Convert the name column to numeric so it will sort properly
        dplyr::mutate(name = as.numeric(name)) %>%

        # Sort by the motif position
        dplyr::arrange(name) %>%

        # Keep only one row with the seed probe z-score
        dplyr::distinct(value, .keep_all = TRUE) %>%

        # Shuffle the z-score column
        dplyr::mutate(value = sample(value))

    # Figure out the new seed probe z-score
    new_seed_probe_zscore <-
        shuffled_zscore_motif %>%

        # Keep only the row with the seed nucleotide in the first position
        dplyr::filter(SNV_nuc == first_seed_nuc & name == 1) %>%

        # Pull the z-score for this row
        dplyr::pull(value)

    # Fill in missing values in the shuffled z-score motif and reformat
    shuffled_zscore_motif <-
        shuffled_zscore_motif %>%

        # Fill in missing SNV_pos_offset/SNV_nuc combos with the seed z-score
        tidyr::complete(
            SNV_nuc,
            name,
            fill = list(value = new_seed_probe_zscore)
        ) %>%

        # Reformat so the columns are positions and the rows are nucleotides
        tidyr::pivot_wider(
            names_from = name,
            values_from = value
        ) %>%

        # Make the nucleotides the row names
        tibble::column_to_rownames("SNV_nuc")

    # Return the shuffled z-score motif
    return(shuffled_zscore_motif)
}


motif_match_pvalue <-
    function(
        corec_motif,
        tf_motif,
        method = c("EUCL", "PCC", "ALLR")
    ) {
        # Make sure the selected comparison method is a valid option
        method <-
            match.arg(method)

        # Make a bunch of shuffled versions of the z-score motif
        shuffled_motifs <-
            replicate(
                1000,
                shuffle_zscore_motif(corec_motif@zscore_motif),
                simplify = FALSE
            )

        # Convert the shuffled z-score motifs into PPMs
        shuffled_ppms <-
            purrr::map(
                shuffled_motifs,
                zscore_to_ppm,
                beta = corec_motif@beta
            )

        # Compare all the shuffled PPMs to the TF motif
        motif_comparisons <-
            universalmotif::compare_motifs(
                c(shuffled_ppms, tf_motif),
                compare.to = length(shuffled_ppms) + 1,
                method = method,
                max.p = Inf,
                max.e = Inf
            ) %>%

            # Don't tell me some comparisons failed I don't care leave me alone
            suppressWarnings() %>%

            # Convert to a data frame
            as.data.frame() %>%

            # Extract just the column of motif comparison scores
            dplyr::pull(score)

        # Calculate the score of the original motif
        original_score <-
            universalmotif::compare_motifs(
                c(corec_motif@ppm, tf_motif),
                compare.to = 2,
                method = method,
                max.p = Inf,
                max.e = Inf
            ) %>%

            # Don't tell me some comparisons failed I don't care leave me alone
            suppressWarnings() %>%

            # Convert to a data frame
            as.data.frame() %>%

            # Extract just the column of motif comparison scores
            dplyr::pull(score)

        # Calculate the one-tailed p-value
        p_val <-
            1 - pnorm(
                original_score,
                mean(motif_comparisons),
                sd(motif_comparisons)
            )

        # Return the p-value
        return(p_val)
    }


identify_motif_match_old <-
    function(
        corec_motif,
        reference_motifs,
        method = c("EUCL", "PCC", "ALLR")
    ) {
        # Make sure the selected comparison method is a valid option
        method <-
            match.arg(method)

        # Compare the corecmotif to the full library of reference motifs
        motif_comparison <-
            universalmotif::compare_motifs(
                c(reference_motifs, corec_motif@ppm),
                compare.to = length(reference_motifs) + 1,
                method = method,
                max.p = Inf,
                max.e = Inf
            ) %>%

            suppressWarnings()

        # Figure out which reference motif has the lowest p-value
        best_match <-
            motif_comparison %>%

            # Convert to a data frame first to remove extra information
            as.data.frame() %>%

            # Take the single row with the smallest p-value
            dplyr::slice_min(Pval, with_ties = FALSE)

        # Update the motif matching slots of the original corecmotif object
        corec_motif@motif_match <-
            reference_motifs[[best_match$target.i]]

        corec_motif@motif_match_score_type <-
            method

        corec_motif@motif_match_score <-
            best_match$score

        corec_motif@motif_match_pvalue <-
            best_match$Pval

        # Return the updated corecmotif
        return(corec_motif)
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

