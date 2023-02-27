#' Import fluorescence data
#'
#' Reads a file containing fluorescence data (in which rows correspond to probes
#' and columns correspond to fluorescence values) and reformats it. The output
#' from this function is suitable for use by [annotate_fluorescence_table].
#'
#' The fluorescence data file supplied to `fluorescence_file` must contain
#' \eqn{n + 3} tab-delimited columns, where n is the number of PBM conditions
#' supplied to `pbm_conditions`. This file must NOT have a header. The columns
#' must be as follows:
#'
#' 1. The probe ID. The orientation tags (i.e., "_o1" or "_o2") will be removed
#' from these IDs.
#' 2. The probe sequence. This column will be removed.
#' 3. The number of PBM conditions profiled on this array. This column will be
#' removed.
#'
#' Columns 4 through \eqn{(n + 3)} must contain the fluorescence values for the
#' PBM conditions in the order they appear in `pbm_conditions`.
#'
#' @param fluorescence_file `character(1)`. The path to the fluorescence data
#'   file to load. See 'Details' for expected columns.
#' @param pbm_conditions `character`. The PBM conditions (e.g., cell type,
#'   treatment, and factor profiled) in the order they appear in
#'   `fluorescence_file`.
#'
#' @return A data frame containing the fluorescence values from the specified
#'   fluorescence data file. The first column will be named "probe_id", and the
#'   remaining columns will be named according to the PBM conditions supplied to
#'   `pbm_conditions`.
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN")
load_fluorescence_data <- function(fluorescence_file, pbm_conditions) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        assertthat::is.string(fluorescence_file) &&
            file.exists(fluorescence_file),
        is.character(pbm_conditions)
    )

    # Load the table of fluorescence values
    fluorescence_table <-
        utils::read.table(
            fluorescence_file,
            header = FALSE,
            sep = "\t",
            strip.white = TRUE,
            stringsAsFactors = FALSE
        )

    # Make sure the list of PBM conditions is the right length
    if (ncol(fluorescence_table) != (length(pbm_conditions) + 3)) {
        stop(
            "pbm_conditions is the wrong length\n",
            "Expected ",
            ncol(fluorescence_table) -3,
            " values but got ",
            length(pbm_conditions),
            call. = FALSE
        )
    }

    # Reformat the fluorescence table
    fluorescence_table <-
        fluorescence_table %>%

        # Rename the columns
        magrittr::set_colnames(
            c("probe_id", "probe_sequence", "n_conditions", pbm_conditions)
        ) %>%

        # Remove the unnecessary columns
        dplyr::select(-probe_sequence, -n_conditions) %>%

        # Remove the orientation tag from the probe IDs
        dplyr::mutate(probe_id = gsub("_o[1-2]", "", probe_id))

    # Return the fluorescence table
    return(fluorescence_table)
}

#' Make a data frame summarizing a list of CoRecMotifs
#'
#' Creates a data frame representation of a list of
#' [CoRecMotifs][CoRecMotif-class].
#'
#' @param corecmotifs `list`. The [CoRecMotifs][CoRecMotif-class] to summarize.
#'
#' @return A data frame with information about a list of
#'   [CoRecMotifs][CoRecMotif-class].
#'
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
#' Updates the `match_cluster` slot of a [CoRecMotif][CoRecMotif-class] based on
#' the provided cluster assignments and the name of the motif in the
#' CoRecMotif's `match_motif` slot.
#'
#' @param corecmotif [CoRecMotif][CoRecMotif-class]. The
#'   [CoRecMotif][CoRecMotif-class] to update.
#' @param cluster_assignments `data.frame` or `NULL`. A table mapping the
#'   reference motifs to motif clusters or NULL to reset the `match_cluster`
#'   slot to NA. See [motif_clusters] for expected columns. (Default: NULL)
#'
#' @return A [CoRecMotif][CoRecMotif-class] with its `match_cluster` slot
#'   updated.
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN")
update_cluster_match <- function(corecmotif, cluster_assignments = NULL) {
    # Clear the cluster match slot if there are no clusters or no motif match
    if (is.null(cluster_assignments) || is.null(match_motif(corecmotif))) {
        corecmotif@match_cluster <- NA_character_

        # Return the updated corecmotif
        return(corecmotif)
    }

    # Clear the cluster match slot if the motif match isn't in the clusters
    if (!(match_altname(corecmotif) %in% cluster_assignments$motif)) {
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
        dplyr::filter(motif == match_altname(corecmotif)) %>%

        # Pull out the cluster name for this motif
        dplyr::pull(cluster)

    # Update the cluster match slot
    corecmotif@match_cluster <- as.character(best_cluster)

    # Return the updated corecmotif
    return(corecmotif)
}

