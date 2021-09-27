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

