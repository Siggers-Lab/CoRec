get_seed_name <- function(corecmotif) {
    return(corecmotif@seed_name)
}

get_pbm_condition <- function(corecmotif) {
    return(corecmotif@pbm_condition)
}

get_array_id <- function(corecmotif) {
    return(corecmotif@array_id)
}

get_beta <- function(corecmotif) {
    return(corecmotif@beta)
}

get_zscore_motif <- function(corecmotif) {
    return(corecmotif@zscore_motif)
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
get_delta_zscore_motif <- function(corecmotif) {
    # Get the z-score motif
    zscore_motif <- get_zscore_motif(corecmotif)

    # Transform the z-scores to reflect their deviation from column-wise median
    delta_zscore_motif <-
        # Subtract the column-wise median from each value in each column
        apply(zscore_motif, 2, function(col) {col - median(col)})

    # Return the delta z-score motif
    return(delta_zscore_motif)
}

get_rolling_ic <- function(corecmotif) {
    return(corecmotif@rolling_ic)
}

get_motif_strength <- function(corecmotif) {
    return(corecmotif@motif_strength)
}

get_seed_sequence <- function(corecmotif) {
    return(corecmotif@seed_sequence)
}

# Seed probe z-score
#
# Finds the fluorescence value z-score of the seed probe of the z-score motif.
#
# @param zscore_motif A data frame representing a z-score motif, where the rows
#   are nucleotides and the columns are positions in the motif.
#
# @return The fluorescence value z-score of the seed probe of the z-score motif.
get_seed_zscore <- function(corecmotif) {
    # Get the z-score motif
    zscore_motif <- get_zscore_motif(corecmotif)

    seed_zscore <- find_seed_zscore(zscore_motif)

    # Return the z-score
    return(seed_zscore)
}

get_ppm <- function(corecmotif) {
    return(corecmotif@ppm)
}

get_motif_name <- function(corecmotif) {
    return(corecmotif@ppm@name)
}

get_match_motif <- function(corecmotif) {
    return(corecmotif@match)
}

get_match_name <- function(corecmotif) {
    return(corecmotif@match@name)
}

get_match_altname <- function(corecmotif) {
    return(corecmotif@match@altname)
}

get_match_pvalue <- function(corecmotif) {
    return(corecmotif@match_pvalue)
}

get_match_cluster <- function(corecmotif) {
    return(corecmotif@match_cluster)
}

#' Make a data frame summarizing a list of CoRecMotifs
#'
#' Creates a data frame representation of a list of
#' [CoRecMotifs][CoRecMotif-class].
#'
#' @param corecmotifs the list of [CoRecMotifs][CoRecMotif-class] to summarize.
#'
#' @return A data frame with information about a list of
#'   [CoRecMotifs][CoRecMotif-class].
#' @export
#'
#' @examples
#' print("FILL THIS IN")
summarize_corecmotifs <- function(corecmotifs) {
    # Convert each corecmotif object into a data frame
    corecmotif_df <-
        lapply(corecmotifs, as.data.frame) %>%

        # Combine all the data frames
        dplyr::bind_rows()

    # Return the dataframe of corecmotif information
    return(corecmotif_df)
}

#' Update the match cluster of a CoRecMotif
#'
#' Update the `match_cluster` slot of a [CoRecMotif][CoRecMotif-class] based on
#' the provided cluster assignments and the name of the motif in the
#' CoRecMotif's `match_motif` slot.
#'
#' @param corecmotif the [CoRecMotif][CoRecMotif-class] to update.
#' @param cluster_assignments a data frame of reference motif names and the
#'   clusters they are assigned to or NULL to reset the `match_cluster` slot to
#'   NA. (Default: NULL)
#'
#' @return A [CoRecMotif][CoRecMotif-class] with its `match_cluster` slot
#'   updated.
#' @export
#'
#' @examples
#' print("FILL THIS IN")
update_cluster_match <- function(corecmotif, cluster_assignments = NULL) {
    # Clear the cluster match slot if there are no clusters or no motif match
    if (is.null(cluster_assignments) || is.null(corecmotif@match_motif)) {
        corecmotif@match_cluster <- NA_character_

        # Return the updated corecmotif
        return(corecmotif)
    }

    # Clear the cluster match slot if the motif match isn't in the clusters
    if (!(corecmotif@match_motif@altname %in% cluster_assignments$motif)) {
        corecmotif@match_cluster <- NA_character_

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
        dplyr::filter(motif == corecmotif@match_motif@altname) %>%

        # Pull out the cluster name for this motif
        dplyr::pull(cluster)

    # Update the cluster match slot
    corecmotif@match_cluster <-
        as.character(best_cluster)

    # Return the updated corecmotif
    return(corecmotif)
}


