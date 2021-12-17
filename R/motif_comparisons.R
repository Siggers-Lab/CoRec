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
        meme_path = "/share/pkg.7/meme/5.3.3/install/bin/",
        min_overlap = 6
    ) {
        # Make sure the selected comparison method is a valid option
        method <-
            match.arg(method)

        # Compare the corecmotif to the full library of reference motifs
        motif_comparison <-
            memes::runTomTom(
                corec_motif@ppm,
                database = reference_motifs_file,
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

        corec_motif@motif_match_score_type <-
            method

        corec_motif@motif_match_score <-
            NA_real_

        corec_motif@motif_match_pvalue <-
            as.numeric(motif_comparison$best_match_pval)

        corec_motif@motif_match_qvalue <-
            as.numeric(motif_comparison$best_match_qval)

        # Return the updated corecmotif
        return(corec_motif)
    }


cluster_motifs <- function(motifs, k = 10, h = NULL) {
    # Compute a distance matrix comparing all the motifs
    compared_motifs <-
        universalmotif::compare_motifs(motifs, method = "EUCL")

    # Convert the distance matrix into a dist object
    compared_motifs <-
        as.dist(compared_motifs)

    # Perform hierarchical clustering on the motifs
    compared_motifs <-
        hclust(compared_motifs)

    # Assign each motif to a cluster by cutting the hierarchical tree
    cluster_assignments <-
        cutree(compared_motifs, k = k, h = h)

    # Make a list of lists of motifs grouped by cluster
    clustered_motifs <-
        lapply(1:max(cluster_assignments), function(cluster) {
            # Find the indices of all the occurrences of this cluster number
            motif_indices <- which(cluster_assignments == cluster)

            # Extract the motifs at those indices
            cluster_motifs <- motifs[motif_indices]
        })

    # Return the list of clusters of motifs
    return(clustered_motifs)
}


assign_to_cluster <- function(corec_motif, clustered_motifs) {
    cluster_scores <-
        lapply(clustered_motifs, function(cluster) {
            # Compare corec_motif to each motif in this cluster
            comparisons <-
                universalmotif::compare_motifs(
                    c(cluster, corec_motif@ppm),
                    compare.to = length(cluster) + 1,
                    method = "EUCL",
                    max.p = Inf,
                    max.e = Inf
                ) %>%

                # Convert output to a data frame
                as.data.frame()

            # Find the average comparison score for this cluster
            return(mean(comparisons$score))
        }) %>%

        # Convert average cluster scores to a vector
        unlist()

    # Find the cluster with the smallest average score
    best_cluster <- which(cluster_scores == min(cluster_scores))

    # Return that cluster number
    return(best_cluster)
}

