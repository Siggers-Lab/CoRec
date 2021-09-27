#' Title
#'
#' Identifies the reference motif that is the best match to the corecmotif and
#' calculates a p-value for that match.
#'
#' @param corec_motif
#' @param reference_motifs
#' @param method
#'
#' @return
#' @export
#'
#' @examples
identify_motif_match <-
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

