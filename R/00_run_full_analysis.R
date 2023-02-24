#' Make a list of CoRecMotifs
#'
#' Creates a list of [CoRecMotifs][CoRecMotif-class] using fluorescence data
#' from a CoRec experiment. This is a convenience function that calls
#' [annotate_fluorescence_table()], [fluorescence_to_zscore_table()], and then
#' [zscore_table_to_corecmotifs()].
#'
#' By default no output files are created. To save output files, you must
#' provide `output_directory`, `output_base_name`, or both. The output files are
#' named according to the following rules.
#'
#' * If `output_directory` is provided but `output_base_name` is not, the files
#'   will be saved in the provided directory using the base name "output".
#' * If `output_base_name` is provided but `output_directory` is not, the files
#'   will be saved in the working directory using the provided base name.
#' * If `array_id` is provided, it will be appended to the output base name.
#' * The individual output file names will be identified by a suffix appended
#'   to the output base name after the array ID (if present).
#' * The output of [annotate_fluorescence_table()] will have the suffix
#'   "fluorescence.tsv".
#' * The output of [fluorescence_to_zscore_table()] will have the suffix
#'   "zscores.tsv"
#' * The output of [zscore_table_to_corecmotifs()] will have the suffix
#'   "all_corecmotifs.rds"
#'
#' @inheritParams annotate_fluorescence_table
#' @inheritParams zscore_table_to_corecmotifs
#' @param fluorescence_file the path to the file containing the fluorescence
#'   data to load. See 'Details' of [annotate_fluorescence_table()] for expected
#'   columns.
#' @param output_directory the directory where output files will be saved. No
#'   output files will be created unless output_directory or output_base_name is
#'   provided. (Default: NULL)
#' @param output_base_name the base name for output files. No output files will
#'   be created unless output_directory or output_base_name is provided.
#'   (Default: NULL)
#'
#' @return A list of [CoRecMotifs][CoRecMotif-class], one for each possible
#'   combination of the probe sets annotated in `annotation` and the PBM
#'   conditions listed in `pbm_conditions`.
#'
#' @seealso [annotate_fluorescence_table()], [fluorescence_to_zscore_table()],
#'   [zscore_table_to_corecmotifs()]
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN!")
make_corecmotifs <-
    function(
        fluorescence_file,
        pbm_conditions,
        annotation,
        array_id = NULL,
        output_directory = NULL,
        output_base_name = NULL
    ) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        assertthat::is.string(fluorescence_file) &&
            file.exists(fluorescence_file),
        is.character(pbm_conditions),
        is.data.frame(annotation),
        assertthat::is.string(output_directory) || is.null(output_directory),
        assertthat::is.string(output_base_name) || is.null(output_base_name),
        assertthat::is.string(array_id) || is.null(array_id)
    )

    # Update the output base name with the output directory and array ID
    output_base_name <-
        update_output_base_name(output_directory, output_base_name, array_id)

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
            annotation = annotation,
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

#' Filter a list of CoRecMotifs and match them to reference motifs
#'
#' Removes [CoRecMotifs][CoRecMotif-class] with low z-scores or motif strengths
#' and CoRecMotifs that do not replicate, then compares the remaining
#' CoRecMotifs to a database of reference motifs to identify the best match.
#' This is a convenience function that calls [filter_corecmotifs()],
#' [check_replicates()], [find_match()], then [filter_corecmotifs()] and
#' [check_replicates()] again, and finally [summarize_corecmotifs()].
#'
#' By default no output files are created. To save output files, you must
#' provide `output_directory`, `output_base_name`, or both. The output files are
#' named according to the following rules.
#'
#' * If `output_directory` is provided but `output_base_name` is not, the files
#'   will be saved in the provided directory using the base name "output".
#' * If `output_base_name` is provided but `output_directory` is not, the files
#'   will be saved in the working directory using the provided base name.
#' * The individual output file names will be identified by a suffix appended
#'   to the output base name.
#' * The output after the first call to [filter_corecmotifs()] and
#'   [check_replicates()] will have the suffix "filtered_corecmotifs.rds".
#' * The output of [find_match()] will have the suffix
#'   "matched_corecmotifs.rds".
#' * The output after the second call to [filter_corecmotifs()] and
#'   [check_replicates()] will have the suffix "significant_corecmotifs.rds".
#' * The output of [summarize_corecmotifs()] will have the suffix
#'   "significant_corecmotifs_summary.tsv".
#'
#' @inheritParams filter_corecmotifs
#' @inheritParams check_replicates
#' @inheritParams find_match
#' @inheritParams make_corecmotifs
#' @param corecmotifs The list of [CoRecMotifs][CoRecMotif-class] to process.
#'
#' @return A filtered list of replicated [CoRecMotifs][CoRecMotif-class] that
#'   match a reference motif.
#'
#' @seealso [filter_corecmotifs()], [check_replicates()], [find_match()],
#'   [summarize_corecmotifs()]
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN")
process_corecmotifs <-
    function(
        corecmotifs,
        reference_motifs_file,
        cluster_assignments = NULL,
        meme_path = NULL,
        rolling_ic = 1,
        motif_strength = 1,
        n_replicates = 2,
        eucl_distance = 0.4,
        min_overlap = 5,
        match_pvalue = 0.05,
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
        update_output_base_name(output_directory, output_base_name)

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
            rolling_ic = rolling_ic,
            motif_strength = motif_strength
        )

    # Filter out CoRecMotifs that don't replicate
    replicated_corecmotifs <-
        check_replicates(
            corecmotifs = filtered_corecmotifs,
            n_replicates = n_replicates,
            eucl_distance = eucl_distance,
            output_file = filtered_output
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
            match_pvalue = match_pvalue
        ) %>%

        # Make sure at least min_n_replicates match a reference motif well
        check_replicates(
            n_replicates = n_replicates,
            eucl_distance = NULL,
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
