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
        reference_motifs_file,
        cluster_assignments = NULL,
        min_overlap = 5,
        method = c(
            "ed",
            "allr",
            "pearson",
            "sandelin",
            "kullback",
            "blic1",
            "blic5",
            "llr1",
            "llr5"
        ),
        meme_path = "/share/pkg.7/meme/5.3.3/install/bin/"
    ) {
        # Make sure the selected comparison method is a valid option
        method <-
            match.arg(method)

        # Compare the corecmotif to the full library of reference motifs
        motif_comparison <-
            memes::runTomTom(
                corec_motif@ppm,
                database = reference_motifs_file,
                thresh = 1,
                evalue = FALSE,
                min_overlap = min_overlap,
                dist = method,
                meme_path = meme_path
            )

        # If no matches were found, return the original corecmotif object
        if (is.na(motif_comparison$best_match_motif)) {
            return(corec_motif)
        }

        # Update the motif matching slots of the original corecmotif object
        corec_motif@motif_match <-
            motif_comparison$best_match_motif[[1]]

        corec_motif@motif_match_method <-
            method

        corec_motif@motif_match_pvalue <-
            as.numeric(motif_comparison$best_match_pval)

        corec_motif@motif_match_qvalue <-
            as.numeric(motif_comparison$best_match_qval)

        # If there are no cluster assignments, return the updated corecmotif
        if (is.null(cluster_assignments)) {
            return(corec_motif)
        }

        # Figure out the cluster the best match motif is in
        best_cluster <-
            cluster_assignments %>%

            # Keep just the row corresponding to the best match motif
            dplyr::filter(
                motif == motif_comparison$best_match_altname
            ) %>%

            # Pull out the cluster name for this motif
            dplyr::pull(cluster)

        # Update the cluster match slot of the original corecmotif object
        corec_motif@motif_cluster_match <-
            as.character(best_cluster)
        
        # Return the updated corecmotif
        return(corec_motif)
    }

