#' Match a list of CoRecMotifs to reference motifs
#'
#' Identifies the reference motif that is the best match to each
#' [CoRecMotif][CoRecMotif-class].
#'
#' The PPM representation of each [CoRecMotif][CoRecMotif-class] is compared to
#' the reference motifs with [memes::runTomTom()] using Euclidean distance
#' ("ed") as the distance metric. The match_motif and match_pvalue slots of the
#' [CoRecMotif][CoRecMotif-class] are filled in based on the best match returned
#' by [memes::runTomTom()]. The match_pvalue slot is corrected for multiple
#' comparisons by multiplying the raw match p-value by the number of motifs in
#' the reference database. The match_cluster slot is filled in based on the
#' user-provided cluster assignments or left empty if no cluster assignments are
#' provided.
#'
#' @inheritParams annotate_fluorescence_table
#' @param corecmotifs `list`. The [CoRecMotifs][CoRecMotif-class] to match to
#'   reference motifs.
#' @param reference_motifs_file `character(1)`. The path to the MEME format file
#'   of reference motifs to match to.
#' @param cluster_assignments `data.frame` or `NULL`. A table mapping the
#'   reference motifs to motif clusters or NULL to skip the cluster assignment
#'   step. See [motif_clusters] for expected columns. (Default: NULL)
#' @param meme_path `character(1)` or `NULL`. The path to "meme/bin/" or NULL to
#'   rely on [memes::runTomTom()] to find it. (Default: NULL)
#' @param min_overlap `integer(1)`. The minimum amount of overlap to require
#'   when comparing a [CoRecMotif][CoRecMotif-class] to a reference motif.
#'   (Default: 5)
#'
#' @return A list of [CoRecMotifs][CoRecMotif-class] with the match_motif,
#'   match_pvalue, and (optionally) match_cluster slots filled in.
#'
#' @seealso [motif_clusters] for a description of the expected columns of
#'   `cluster_assignments`
#' @export
#'
#' @examples
#' print("FILL THIS IN")
find_match <-
    function(
        corecmotifs,
        reference_motifs_file,
        cluster_assignments = NULL,
        meme_path = NULL,
        min_overlap = 5,
        output_file = NULL
    ) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        assertthat::is.string(reference_motifs_file) &&
            file.exists(reference_motifs_file),
        is.data.frame(cluster_assignments) || is.null(cluster_assignments),
        assertthat::is.string(meme_path) || is.null(meme_path),
        assertthat::is.count(min_overlap),
        assertthat::is.string(output_file) || is.null(output_file)
    )

    # Make sure corecmotifs is a valid list of CoRecMotifs
    corecmotifs <- check_corecmotif_list(corecmotifs)

    # Make sure cluster_assignments has the expected columns and remove extras
    if (!is.null(cluster_assignments)) {
        cluster_assignments <-
            check_colnames(cluster_assignments, c("motif", "cluster"))
    }

    # Make sure the provided meme_path is valid
    if (!memes::meme_is_installed(meme_path)) {
        stop(
            "Invalid meme_path: could not find MEME installation.",
            call. = FALSE
        )
    }

    # Get the PPMs of the CoRecMotifs to compare
    motifs <- lapply(corecmotifs, get_motif)

    # Compare the PPMs to the full library of reference motifs
    motif_comparison <-
        memes::runTomTom(
            motifs,
            database = reference_motifs_file,
            thresh = 1,
            evalue = FALSE,
            min_overlap = min_overlap,
            dist = "ed",
            meme_path = meme_path
        )

    # Update the match slots of the CoRecMotifs
    corecmotifs <- lapply(1:length(corecmotifs), function(index) {
        # Set the match motif
        corecmotifs[[index]]@match_motif <-
            motif_comparison[[index]]$best_match_motif[[1]]

        # Save the match orientation in the extrainfo slot
        corecmotifs[[index]]@match_motif["extrainfo"] <-
            motif_comparison[[index]]$best_match_strand[[1]]

        # Set the match p-value (actually an E-value, but whatever)
        corecmotifs[[index]]@match_pvalue <-
            as.numeric(motif_comparison[[index]]$best_match_eval)

        # Figure out the match cluster
        corecmotifs[[index]] <-
            update_cluster_match(corecmotifs[[index]], cluster_assignments)

        # Return the updated CoRecMotif
        return(corecmotifs[[index]])
    })

    # Try to save the matched CoRecMotifs as an RDS file if necessary
    try_catch_save_output(corecmotifs, output_file, "rds")

    # Return the list of updated CoRecMotifs
    return(corecmotifs)
}

