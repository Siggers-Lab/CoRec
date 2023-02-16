# Import the universalmotif class to use it in a class union
#' @importClassesFrom universalmotif universalmotif

# Define a class union that can contain either a universalmotif object or NULL
setClassUnion("universalmotif_or_NULL", c("universalmotif", "NULL"))


#' Title
#'
#' @slot seed_name character.
#' @slot pbm_condition character.
#' @slot beta numeric.
#' @slot zscore_motif data.frame.
#' @slot delta_zscore_motif data.frame.
#' @slot ppm universalmotif.
#' @slot seed_zscore numeric.
#' @slot rolling_ic numeric.
#' @slot motif_strength numeric.
#' @slot seed_probe_sequence character.
#' @slot motif_match universalmotif_or_NULL.
#' @slot motif_match_method character.
#' @slot motif_match_pvalue numeric.
#' @slot motif_match_qvalue numeric.
#' @slot motif_cluster_match character.
#'
#' @return
#' @export
#'
#' @examples
setClass(
    # Name the class CoRecMotif
    "corecmotif",

    # Define the names and types of the slots the class should have
    slots = list(
        seed_name = "character",
        pbm_condition = "character",
        beta = "numeric",
        zscore_motif = "data.frame",
        delta_zscore_motif = "data.frame",
        ppm = "universalmotif",
        seed_zscore = "numeric",
        rolling_ic = "numeric",
        motif_strength = "numeric",
        seed_probe_sequence = "character",
        motif_match = "universalmotif_or_NULL",
        motif_match_method = "character",
        motif_match_pvalue = "numeric",
        motif_match_qvalue = "numeric",
        motif_cluster_match = "character"
    ),

    # Provide a default example object
    prototype = list(
        seed_name = NA_character_,
        pbm_condition = NA_character_,
        beta = NA_real_,
        zscore_motif = data.frame(),
        delta_zscore_motif = data.frame(),
        ppm = universalmotif::create_motif("ACGT"),
        seed_zscore = NA_real_,
        rolling_ic = NA_real_,
        motif_strength = NA_real_,
        seed_probe_sequence = NA_character_,
        motif_match = NULL,
        motif_match_method = NA_character_,
        motif_match_pvalue = NA_real_,
        motif_match_qvalue = NA_real_,
        motif_cluster_match = NA_character_
    )
) %>%

    # Do not print the message in universalmotif's constructor
    suppressMessages()



#' Title
#'
#' @param seed_name
#' @param pbm_condition
#' @param beta
#' @param zscore_motif
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
corecmotif <-
    function(
        seed_name,
        pbm_condition,
        zscore_motif,
        ic_window_width = 5,
        top_n_percent = 15,
        ...
    ) {
        delta_zscore_motif <-
            zscore_to_delta_zscore(zscore_motif)

        motif_name <-
            paste(seed_name, pbm_condition, sep = "_")

        beta <-
            calculate_beta(zscore_motif)

        ppm <-
            zscore_to_ppm(zscore_motif, beta, motif_name)

        seed_zscore <-
            find_seed_zscore(zscore_motif)

        rolling_ic <-
            calculate_rolling_ic(ppm, ic_window_width)

        motif_strength <-
            calculate_strength(zscore_motif, top_n_percent)

        new(
            "corecmotif",
            seed_name = seed_name,
            pbm_condition = pbm_condition,
            beta = beta,
            zscore_motif = zscore_motif,
            delta_zscore_motif = delta_zscore_motif,
            ppm = ppm,
            seed_zscore = seed_zscore,
            rolling_ic = rolling_ic,
            motif_strength = motif_strength,
            ...
        )
    }


# Calculate beta
#
# Calculates the beta parameter to use when converting z-score motifs to PPMs.
#
# @param zscore_motif A data frame representing a z-score motif, where the rows
#   are nucleotides and the columns are positions in the motif.
#
# @return The value of beta for the given PBM condition.
calculate_beta <- function(zscore_motif) {
    # Find the seed probe z-score
    seed_probe <- find_seed_zscore(zscore_motif)

    # Calculate beta
    beta <- 4 - (0.5 * seed_probe)

    # Restrict beta to a range of 1 to 4
    beta <- max(min(4, beta), 1)

    # Return beta
    return(beta)

    # Return beta
    return(beta)
}


# Convert a z-score motif to a delta z-score motif
#
# Transforms the given z-score motif by subtracting the median z-score at each
# position from the z-score for each probe at that position.
#
# @param zscore_motif A data frame representing a z-score motif, where the rows
#   are nucleotides and the columns are positions in the motif.
#
# @return A data frame representing the delta z-score motif corresponding to the
#   provided z-score motif, where the rows are nucleotides and the columns are
#   positions in the motif.
zscore_to_delta_zscore <- function(zscore_motif) {
    # Transform the z-scores to reflect their deviation from column-wise median
    delta_zscore_motif <-
        zscore_motif %>%

        # Subtract the column-wise median from each value in each column
        dplyr::mutate_all(list(~ . - median(.)))

    # Return the delta z-score motif
    return(delta_zscore_motif)
}


# Convert a z-score motif to a PPM
#
# Transforms the given z-score motif into a Position Probability Matrix (PPM).
#
# @param zscore_motif A data frame representing a z-score motif, where the rows
#   are nucleotides and the columns are positions in the motif.
# @param beta A normalization factor that scales with the maximum z-score
#   observed for a given PBM condition.
#
# @return A universalmotif object representing the PPM corresponding to the
#   given z-score motif.
zscore_to_ppm <- function(zscore_motif, beta, name = "motif") {
    # Transform the z-scores using the beta parameter
    ppm <-
        # Multiply each z-score by beta and then take the exponential
        exp(beta * zscore_motif) %>%

        # Normalize by dividing each value in each column by the column-wise sum
        dplyr::mutate_all(list(~ . / sum(.))) %>%

        # Convert to a numeric matrix
        as.matrix() %>%

        # Convert to a universalmotif object
        universalmotif::create_motif(name = name)

    # Return the PPM
    return(ppm)
}


calculate_strength <- function(zscore_motif, top_n_percent = 15) {
    # Get a sorted list of all the probe z-scores
    z_scores <-
        sort(unique(unlist(zscore_motif)), decreasing = TRUE)

    # Figure out how many probes to average (rounded up to the nearest integer)
    num_probes <- ceiling(length(z_scores) * (top_n_percent / 100))

    # Find the median z-score of the highest num_probes probes
    median_zscore <- median(z_scores[1:num_probes])

    names(median_zscore) <- top_n_percent

    # Return the average z-score
    return(median_zscore)
}


#' Update the motif cluster match of a corecmotif
#'
#' @param corecmotif
#' @param cluster_assignments
#'
#' @return
#' @export
#'
#' @examples
update_cluster_match <- function(corecmotif, cluster_assignments = NULL) {
    # Clear the cluster match slot if there are no clusters or no motif match
    if (is.null(cluster_assignments) | is.null(corecmotif@motif_match)) {
        corecmotif@motif_cluster_match <- NA_character_

        # Return the updated corecmotif
        return(corecmotif)
    }

    # Clear the cluster match slot if the motif match isn't in the clusters
    if (!(corecmotif@motif_match@altname %in% cluster_assignments$motif)) {
        corecmotif@motif_cluster_match <- NA_character_

        # Print a warning message
        warning(
            "Motif match altname not in cluster assignments table; ",
            "setting cluster match to NA"
        )

        # Return the updated corecmotif
        return(corecmotif)
    }

    # Figure out the cluster the best match motif is in
    best_cluster <-
        cluster_assignments %>%

        # Keep just the row corresponding to the best match motif
        dplyr::filter(
            motif == corecmotif@motif_match@altname
        ) %>%

        # Pull out the cluster name for this motif
        dplyr::pull(cluster)

    # Update the cluster match slot
    corecmotif@motif_cluster_match <-
        as.character(best_cluster)

    # Return the updated corecmotif
    return(corecmotif)
}


# Seed probe z-score
#
# Finds the fluorescence value z-score of the seed probe of the z-score motif.
#
# @param zscore_motif A data frame representing a z-score motif, where the rows
#   are nucleotides and the columns are positions in the motif.
#
# @return The fluorescence value z-score of the seed probe of the z-score motif.
find_seed_zscore <- function(zscore_motif) {
    seed_probe_zscore <-
        # The seed probe z-score shows up at every position of the z-score motif
        which(table(zscore_motif) >= ncol(zscore_motif)) %>%

        # The names of the frequency table are the z-scores
        names() %>%

        # Convert to numeric
        as.numeric()

    # Return the z-score
    return(seed_probe_zscore)
}


# Title
#
# @param ppm
# @param width
#
# @return
calculate_rolling_ic <- function(ppm, width = 5) {
    # Convert the PPM to an information content matrix
    icm <-
        universalmotif::convert_type(ppm, type = "ICM")

    #
    max_sliding_window_ic <-
        icm["motif"] %>%
        colSums() %>%
        zoo::rollmean(width) %>%
        max()

    names(max_sliding_window_ic) <- width

    return(max_sliding_window_ic)
}
