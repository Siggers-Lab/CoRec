#' Make a list of \linkS4class{CoRecMotif} objects
#'
#' Creates a list of \linkS4class{CoRecMotif} objects using fluorescence data
#' from a CoRec experiment.
#'
#' This is a convenience function that calls
#' \code{\link{annotate_fluorescence_table}},
#' \code{\link{fluorescence_to_zscore_table}}, and
#' \code{\link{zscore_table_to_corecmotifs}} and optionally saves the outputs
#' using file names generated from user input. By default no output files are
#' created. To save output files, you must provide \code{output_directory},
#' \code{output_base_name}, or both. If \code{output_directory} is provided but
#' \code{output_base_name} is not, the files will be saved in the provided
#' directory using the base name "output". If \code{output_base_name} is
#' provided but \code{output_directory} is not, the files will be saved in the
#' current working directory using the provided base name. If \code{array_id} is
#' provided, it will be appended to \code{output_base_name}. The outputs of
#' \code{\link{annotate_fluorescence_table}},
#' \code{\link{fluorescence_to_zscore_table}}, and
#' \code{\link{zscore_table_to_corecmotifs}} will be saved with the suffixes
#' "_fluorescence.tsv", "_zscores.tsv", and "_all_corecmotifs.rds" respectively.
#'
#' @param fluorescence_file the path to the file containing the fluorescence
#'   data to load. See 'Details' of \code{\link{annotate_fluorescence_table}}
#'   for expected columns.
#' @param pbm_conditions a character vector specifying the PBM conditions (e.g.,
#'   cell type, treatment, and factor profiled) in the order they appear in
#'   \code{fluorescence_file}.
#' @param annotation_file the path to the file containing the probe annotations
#'   to use. See 'Details' of \code{\link{annotate_fluorescence_table}} for
#'   expected columns.
#' @param output_directory the directory where output files will be saved. No
#'   output files will be created unless output_directory or output_base_name is
#'   provided.
#' @param output_base_name the base name for output files. No output files will
#'   be created unless output_directory or output_base_name is provided.
#' @param array_id an optional (but recommended) tag specifying the particular
#'   array/experiment the fluorescence data is from.
#'
#' @return A list of \linkS4class{CoRecMotif} objects, one for each possible
#'   combination of the probe sets annotated in \code{annotation_file} and the
#'   PBM conditions listed in \code{pbm_conditions}.
#'
#' @export
#'
#' @seealso \link{annotate_fluorescence_table},
#'   \link{fluorescence_to_zscore_table}, \link{zscore_table_to_corecmotifs}
#' @examples
#' corecmotifs <-
#'     make_corecmotifs(
#'         fluorescence_file =
#'             "example_data/hTF_v1_example_fluorescence_rep1.dat",
#'         pbm_conditions = c(
#'             "UT_SUDHL4_SWISNF_mix",
#'             "UT_SUDHL4_HDAC1_mix",
#'             "UT_SUDHL4_PRMT5",
#'             "UT_SUDHL4_JMJD2A"
#'         ),
#'         annotation_file = "example_data/hTF_v1_example_annotation.tsv"
#'         array_id = "v1_a6_run1"
#'     )
make_corecmotifs <-
    function(
        fluorescence_file,
        pbm_conditions,
        annotation_file,
        array_id = NULL,
        output_directory = NULL,
        output_base_name = NULL
    ) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        assertthat::is.string(fluorescence_file) &&
            file.exists(fluorescence_file),
        is.character(pbm_conditions),
        assertthat::is.string(annotation_file) &&
            file.exists(annotation_file),
        assertthat::is.string(output_directory) || is.null(output_directory),
        assertthat::is.string(output_base_name) || is.null(output_base_name),
        assertthat::is.string(array_id) || is.null(array_id)
    )

    # Update the output base name with the output directory and array ID
    output_base_name <-
        update_output_base_name(output_base_name, output_directory, array_id)

    # Do not save any output files by default
    fluorescence_output <- NULL
    zscore_output <- NULL
    corec_output <- NULL

    # Create names for the output files if necessary
    if (!is.null(output_base_name)) {
        fluorescence_output <- paste0(output_base_name, "_fluorescence.tsv")
        zscore_output <- paste0(output_base_name, "_zscores.tsv")
        corec_output <- paste0(output_base_name, "_all_corecmotifs.rds")
    }

    # Load and annotate the table of fluorescence values
    fluorescence_table <-
        annotate_fluorescence_table(
            fluorescence_file = fluorescence_file,
            pbm_conditions = pbm_conditions,
            annotation_file = annotation_file,
            output_file = fluorescence_output
        )

    # Convert the fluorescence values into condition-wise z-scores
    zscore_table <-
        fluorescence_to_zscore_table(
            fluorescence_table = fluorescence_table,
            fluorescence_columns = pbm_conditions,
            output_file = zscore_output
        )

    # Make CoRecMotifs for all the seed/condition combos
    corecmotifs <-
        zscore_table_to_corecmotifs(
            zscore_table = zscore_table,
            zscore_columns = pbm_conditions,
            array_id = array_id,
            output_file = corec_output
        )

    # Return the list of CoRecMotifs
    return(corecmotifs)
}

pipeline_part_2 <-
    function(
        corecmotifs,
        reference_motifs_file,
        min_rolling_ic = 1,
        min_motif_strength = 1,
        min_n_replicates = 2,
        max_eucl_distance = 0.4,
        min_overlap = 5,
        cluster_assignments = NULL,
        max_match_pvalue = 0.05,
        meme_path = "/share/pkg.7/meme/5.3.3/install/bin/"
        output_directory = NULL,
        output_base_name = NULL
    ) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        assertthat::is.string(reference_motifs_file) &&
            file.exists(reference_motifs_file),
        assertthat::is.string(output_directory) || is.null(output_directory),
        assertthat::is.string(output_base_name) || is.null(output_base_name),
        is.data.frame(cluster_assignments) || is.null(cluster_assignments),
        assertthat::is.string(meme_path) || is.null(meme_path),
        assertthat::is.number(rolling_ic) || is.null(rolling_ic),
        assertthat::is.number(motif_strength) || is.null(motif_strength),
        assertthat::is.count(n_replicates),
        assertthat::is.number(eucl_distance) || is.null(eucl_distance),
        assertthat::is.count(min_overlap),
        assertthat::is.number(match_pvalue) || is.null(match_pvalue)
    )

    # Update the output base name with the output directory
    output_base_name <-
        update_output_base_name(output_base_name, output_directory)

    # Do not save any output files by default
    filtered_output <- NULL
    matched_output <- NULL
    final_output <- NULL
    summary_output <- NULL

    # Create names for the output files if necessary
    if (!is.null(output_base_name)) {
        filtered_output <-
            paste0(output_base_name, "_filtered_corecmotifs.rds")
        matched_output <-
            paste0(output_base_name, "_matched_corecmotifs.rds")
        final_output <-
            paste0(output_base_name, "_significant_corecmotifs.rds")
        summary_output <-
            paste0(output_base_name, "_significant_corecmotifs_summary.tsv")
    }

    # Filter out CoRecMotifs with low motif strength and/or rolling IC scores
    filtered_corecmotifs <-
        filter_corecmotifs(
            corecmotifs,
            rolling_ic = min_rolling_ic,
            motif_strength = min_motif_strength
        )

    # Filter out CoRecMotifs that don't replicate
    replicated_corecmotifs <-
        match_replicates(
            corecmotifs = filtered_corecmotifs,
            min_n_replicates = min_n_replicates,
            max_eucl_distance = max_eucl_distance
        )

    # Find the best matching reference motif for each CoRecMotif
    matched_corecmotifs <-
        find_match(
            corecmotifs = replicated_corecmotifs,
            reference_motifs_file = reference_motifs_file,
            cluster_assignments = cluster_assignments,
            min_overlap = min_overlap,
            meme_path = meme_path,
            output_file = matched_output
        )

    # Filter out CoRecMotifs that don't match any reference motifs well
    final_corecmotifs <-
        filter_corecmotifs(
            matched_corecmotifs,
            match_pvalue = max_match_pvalue
        ) %>%

        # Make sure at least min_n_replicates match a reference motif well
        match_replicates(
            min_n_replicates = min_n_replicates,
            max_eucl_distance = NULL,
            output_file = final_output
        )

    # Make and save a summary table if necessary
    if (!is.null(summary_output)) {
        corecmotif_summary <-
            summarize_corecmotifs(final_corecmotifs)

        try_catch_save_output(corecmotif_summary, summary_output, "tsv")
    }

    # Return the list of CoRecMotifs
    return(final_corecmotifs)
}

