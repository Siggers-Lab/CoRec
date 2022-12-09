#' Extract corecmotif objects from z-score matrix
#'
#' Extracts a corecmotif object from the input z-score matrix for every possible
#' combination of the given seed names and PBM conditions.
#'
#' @rdname extract_corecmotifs
#'
#' @param zscore_matrix A data frame of PBM condition-wise fluorescence value
#'   z-scores, such as output by fluorescence_to_zscore_matrix().
#' @param seed_names,seed_names A character vector of the names of the seeds for
#'   which to extract corecmotif objects.
#' @param pbm_condition,pbm_conditions A character vector of the column names of
#'   the PBM conditions for which to extract corecmotif objects.
#'
#' @return A list of corecmotif objects, one for each possible combination of
#'   the seed names and PBM conditions provided to \code{seed_names} and
#'   \code{pbm_conditions} respectively.
#' @export
#'
#' @examples
#' corec_motifs <-
#'     zscore_matrix_to_motifs(
#'         zscore_matrix,
#'         seed_names,
#'         pbm_conditions_br
#'     )
zscore_matrix_to_motifs <-
    function(
        zscore_matrix,
        seed_names,
        pbm_conditions
    ) {
        # Make a data frame of the seed names, PBM conditions, and motifs
        zscore_motif_table <-
            # Get all possible combinations of seed names and PBM conditions
            tidyr::crossing(seed_names, pbm_conditions) %>%

            # Calculate the z-score motif for each seed/condition combination
            dplyr::mutate(
                zscore_motif = purrr::map2(
                    seed_names,
                    pbm_conditions,
                    make_zscore_motif,
                    zscore_matrix = zscore_matrix
                )
            )

        # Make a table of the seed names and sequences
        seeds <-
            zscore_matrix %>%

            # Keep only the seed probe row for each TF probe set
            dplyr::filter(SNV_pos_offset == 0) %>%

            # Select just the seed name and probe sequence columns
            dplyr::select(seed_names, probe_seq)

        # Merge the seed sequences into the motif table
        zscore_motif_table <-
            zscore_motif_table %>%

            dplyr::left_join(seeds, by = "seed_names")

        # Convert the data frame into a list of corecmotif objects
        corec_motifs <-
            purrr::pmap(
                list(
                    seed_name = zscore_motif_table$seed_names,
                    pbm_condition = zscore_motif_table$pbm_conditions,
                    zscore_motif = zscore_motif_table$zscore_motif,
                    seed_probe_sequence = zscore_motif_table$probe_seq
                ),
                corecmotif
            )

        # Return the list of corecmotif objects
        return(corec_motifs)
}


#' Extract corecmotif objects from z-score matrix

#' @rdname extract_corecmotifs
#'
#' @param zscore_matrix A data frame of PBM condition-wise fluorescence value
#'   z-scores, such as output by fluorescence_to_zscore_matrix().
#' @param seed_name The name of the seed for which to extract a corecmotif
#'   object.
#' @param pbm_condition The column name of the PBM condition for which to
#'   extract a corecmotif object.
#'
#' @return A data frame representing the z-score motif for the seed name and PBM
#'   condition provided to \code{seed_name} and \code{pbm_condition}
#'   respectively, where the rows are nucleotides and the columns are positions
#'   in the motif.
#' @export
#'
#' @examples
#' zscore_motif <-
#'     make_zscore_motif(
#'         zscore_matrix,
#'         "MA0079.3_SP1",
#'         "v1_a11_run1_br_UT_SUDHL4_SMARCA4MIX"
#'     )
make_zscore_motif <- function(zscore_matrix, seed_name, pbm_condition) {
    # Get the z-score of the seed probe for this seed_name/pbm_condition combo
    seed_zscore <-
        zscore_matrix[
            which(
                zscore_matrix$SNV_pos_offset == 0 &
                    zscore_matrix$seed_names == seed_name),
            pbm_condition
            ]

    # Fill in a motif data frame with the seed and SV probe z-scores
    zscore_motif <-
        zscore_matrix %>%

        # Keep only the SV probe rows from the relevant probe set
        dplyr::filter(seed_names == seed_name & SNV_pos_offset > 0) %>%

        # Create a column replacing NA z-scores in the relevant condition
        dplyr::mutate(
            zscore = tidyr::replace_na(!!as.symbol(pbm_condition), 0.0000001)
        ) %>%

        # Keep only the nucleotide, position, and z-score columns
        dplyr::select(SNV_pos_offset, SNV_nuc, zscore) %>%

        # Fill in missing SNV_pos_offset/SNV_nuc combos with the seed z-score
        tidyr::complete(
            SNV_pos_offset,
            SNV_nuc,
            fill = list(zscore = seed_zscore)
        ) %>%

        # Reformat so the columns are positions and the rows are nucleotides
        tidyr::pivot_wider(
            names_from = SNV_pos_offset,
            values_from = zscore
        ) %>%

        # Make the nucleotides the row names
        tibble::column_to_rownames("SNV_nuc")

    # Return the z-score motif
    return(zscore_motif)
}


# Calculate beta
#
# Calculates the beta parameter to use when converting z-score motifs to PPMs.
#
# @param zscore_motif A data frame representing a z-score motif, where the rows
#   are nucleotides and the columns are positions in the motif.
#
# @return The value of beta for the given PBM condition.
calculate_beta <- function(zscore_motif, method = c("linear", "old")) {
    # Make sure the selected method for calculating beta is a valid option
    method <-
        match.arg(method)

    # Calculate beta using the selected method
    beta <-
        switch(
            method,

            # Old: 15 / the maximum z-score from any probe in this probe set
            "old" = 15 / max(zscore_motif),

            # Linear: between 1 and 4 based on seed probe z-score
            "linear" = {
                # Find the seed probe z-score
                seed_probe <- find_seed_zscore(zscore_motif)

                # Calculate beta
                beta <- 4 - (0.5 * seed_probe)

                # Restrict beta to a range of 1 to 4
                beta <- max(min(4, beta), 1)

                # Return beta
                return(beta)
            }
        )

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

