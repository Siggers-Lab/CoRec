#' Match a list of \linkS4class{CoRecMotif} objects to reference motifs
#'
#' Identifies the reference motif that is the best match to each
#' \linkS4class{CoRecMotif}.
#'
#' The PPM representation of each CoRecMotif is compared to the reference motifs
#' with \link{\code{memes::runTomTom()}} using Euclidean distance ("ed") as the
#' distance metric. The match_motif and match_pvalue slots are filled in based
#' on the best match returned by \link{\code{memes::runTomTom()}}. The
#' match_pvalue slot is corrected for multiple comparisons by multiplying the
#' raw match p-value by the number of motifs in the reference database. The
#' match_cluster slot is filled in based on the user-provided cluster
#' assignments or left empty if no cluster assignments are provided.
#'
#' @param corecmotifs the list of CoRecMotif objects to match to reference
#'   motifs.
#' @param reference_motifs_file the path to a MEME format file of reference
#'   motifs to match to.
#' @param cluster_assignments a data frame mapping the reference motifs to motif
#'   clusters or NULL to skip the cluster assignment step.
#' @param min_overlap a single positive integer specifying the minimum amount of
#'   overlap to require when comparing a CoRecMotif to a reference motif.
#' @param meme_path the path to "meme/bin/" or NULL to rely on
#'   \link{\code{memes::runTomTom()}} to find it.
#'
#' @return A list of \linkS4class{CoRecMotif} objects with the match_motif,
#'   match_pvalue, and (optionally) match_cluster slots filled in.
#'
#' @export
#'
#' @examples
#' # Load example CoRecMotifs
#' corecmotifs <-
#'     readRDS(
#'         "example_data/output/example_rep1_v1_a11_run1_all_corecmotifs.rds"
#'     )
#'
#' # Compare the CoRecMotifs to the reference motifs
#' matched_corecmotifs <-
#'     find_match(
#'         corecmotifs,
#'         "example_data/Homo_sapiens_JASPAR2022_CORE_filtered.meme"
#'     )
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

    # Make sure corecmotifs is a list
    if (!is.list(corecmotifs)) {
        corecmotifs <- list(corecmotifs)
    }

    # Get the PPMs of the CoRecMotifs to compare
    ppms <- lapply(corecmotifs, get_ppm)

    # Compare the PPMs to the full library of reference motifs
    motif_comparison <-
        memes::runTomTom(
            ppms,
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

