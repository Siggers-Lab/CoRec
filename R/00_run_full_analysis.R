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
#' @param output_directory `character(1)` or `NULL`. The path to the directory
#'   where output files will be saved. No output files will be created unless
#'   `output_directory` or `output_base_name` is provided. (Default: NULL)
#' @param output_base_name `character(1)` or `NULL`. The base name for output
#'   files. No output files will be created unless `output_directory` or
#'   `output_base_name ` is provided. (Default: NULL)
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
        fluorescence_table,
        fluorescence_columns,
        annotation,
        array_id = NULL,
        output_directory = NULL,
        output_base_name = NULL
    ) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        is.data.frame(fluorescence_table),
        is.character(fluorescence_columns),
        is.data.frame(annotation),
        assertthat::is.string(array_id) || is.null(array_id),
        assertthat::is.string(output_directory) || is.null(output_directory),
        assertthat::is.string(output_base_name) || is.null(output_base_name)
    )

    # If no array ID is given, generate a random ID
    # CoRecMotifs have to have different names or check_replicates() won't work
    if (is.null(array_id)) {
        array_id <- create_array_id()
    }

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
    annotated_fluorescence_table <-
        annotate_fluorescence_table(
            fluorescence_table = fluorescence_table,
            fluorescence_columns = fluorescence_columns,
            annotation = annotation,
            output_file = fluorescence_output
        )

    # Convert the fluorescence values into condition-wise z-scores
    zscore_table <-
        fluorescence_to_zscore_table(
            fluorescence_table = annotated_fluorescence_table,
            fluorescence_columns = fluorescence_columns,
            output_file = zscore_output
        )

    # Make CoRecMotifs for all the seed/condition combos
    corecmotifs <-
        zscore_table_to_corecmotifs(
            zscore_table = zscore_table,
            zscore_columns = fluorescence_columns,
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
#' [check_replicates()] again, and finally [summarize_corecmotifs()] and
#' optionally [compare_conditions()].
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
#' * The output of [summarize_corecmotifs()] with `by_cluster = TRUE` will have
#'   the suffix "significant_corecmotifs_summary_by_cluster.tsv".
#' * The output of [compare_conditions()] will have the suffix
#'   "condition_comparisons.tsv".
#'
#' @inheritParams check_replicates
#' @inheritParams find_match
#' @inheritParams make_corecmotifs
#' @param corecmotifs `list`. The [CoRecMotifs][CoRecMotif-class] to process.
#' @param motif_strength `numeric(1)` or `NULL`. The minimum motif strength to
#'   keep or NULL not to filter by motif strength. (Default: 1)
#' @param rolling_ic `numeric(1)` or `NULL`. The minimum rolling IC to keep or
#'   NULL not to filter by rolling IC. (Default: 1)
#' @param match_pvalue `numeric(1)` or `NULL`. The maximum match p-value to keep
#'   or NULL not to filter by match p-value. (Default: 0.01)
#' @param pbm_condition_groups `list(character)` or `NULL`. The names of the PBM
#'   conditions to compare. Each element of the list should contain a group of
#'   PBM conditions to compare to each other. If the list elements are named,
#'   the names will be passed to the `pbm_conditions_group` parameter of
#'   [compare_conditions()]. (Default: NULL)
#'
#' @return A filtered list of replicated [CoRecMotifs][CoRecMotif-class] that
#'   match a reference motif.
#'
#' @seealso [filter_corecmotifs()], [check_replicates()], [find_match()],
#'   [summarize_corecmotifs()], [compare_conditions()]
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
        motif_strength = 1,
        rolling_ic = 1,
        n_replicates = 2,
        eucl_distance = 0.4,
        min_overlap = 5,
        match_pvalue = 0.01,
        pbm_condition_groups = NULL,
        output_directory = NULL,
        output_base_name = NULL
    ) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        assertthat::is.string(reference_motifs_file) &&
            file.exists(reference_motifs_file),
        is.data.frame(cluster_assignments) || is.null(cluster_assignments),
        assertthat::is.string(meme_path) || is.null(meme_path),
        assertthat::is.number(motif_strength) || is.null(motif_strength),
        assertthat::is.number(rolling_ic) || is.null(rolling_ic),
        assertthat::is.count(n_replicates),
        assertthat::is.number(eucl_distance) || is.null(eucl_distance),
        assertthat::is.count(min_overlap),
        assertthat::is.number(match_pvalue) || is.null(match_pvalue),
        is.null(pbm_condition_groups) || (
            is.list(pbm_condition_groups) &&
                all(vapply(pbm_condition_groups, is.character, logical(1)))
            ),
        assertthat::is.string(output_directory) || is.null(output_directory),
        assertthat::is.string(output_base_name) || is.null(output_base_name)
    )

    # Update the output base name with the output directory
    output_base_name <-
        update_output_base_name(output_directory, output_base_name)

    # Do not save any output files by default
    filtered_output <- NULL
    matched_output <- NULL
    final_output <- NULL
    summary_output <- NULL
    cluster_output <- NULL
    comparison_output <- NULL

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
        cluster_output <-
            paste0(
                output_base_name,
                "_significant_corecmotifs_summary_by_cluster.tsv"
            )
        comparison_output <-
            paste0(output_base_name, "_condition_comparisons.tsv")
    }

    # Filter out CoRecMotifs with low motif strength and/or rolling IC scores
    filtered_corecmotifs <-
        filter_corecmotifs(
            corecmotifs,
            motif_strength = motif_strength,
            rolling_ic = rolling_ic,
            check_corecmotifs = TRUE
        )

    # Filter out CoRecMotifs that don't replicate
    replicated_corecmotifs <-
        check_replicates(
            corecmotifs = filtered_corecmotifs,
            n_replicates = n_replicates,
            eucl_distance = eucl_distance,
            output_file = filtered_output,
            check_corecmotifs = FALSE
        )

    # Find the best matching reference motif for each CoRecMotif
    matched_corecmotifs <-
        find_match(
            corecmotifs = replicated_corecmotifs,
            reference_motifs_file = reference_motifs_file,
            cluster_assignments = cluster_assignments,
            min_overlap = min_overlap,
            meme_path = meme_path,
            output_file = matched_output,
            check_corecmotifs = FALSE
        )

    # Filter out CoRecMotifs that don't match any reference motifs well
    final_corecmotifs <-
        filter_corecmotifs(
            matched_corecmotifs,
            match_pvalue = match_pvalue,
            check_corecmotifs = FALSE
        ) %>%

        # Make sure at least min_n_replicates match a reference motif well
        check_replicates(
            n_replicates = n_replicates,
            eucl_distance = NULL,
            output_file = final_output,
            check_corecmotifs = FALSE
        )

    # Make and save summary tables if necessary
    if (!is.null(summary_output)) {
        # Summarize at the level of individual motifs
        corecmotif_summary <-
            summarize_corecmotifs(final_corecmotifs, check_corecmotifs = FALSE)

        try_catch_save_output(corecmotif_summary, summary_output, "tsv")

        # Summarize at the level of clusters
        cluster_summary <-
            summarize_corecmotifs(
                final_corecmotifs,
                by_cluster = TRUE,
                check_corecmotifs = FALSE
            )

        try_catch_save_output(cluster_summary, cluster_output, "tsv")

        # Compare across conditions if necessary
        if (!is.null(pbm_condition_groups)) {
            condition_comparison_df <-
                lapply(1:length(pbm_condition_groups), function(i) {
                    compare_conditions(
                        matched_corecmotifs,
                        pbm_conditions = pbm_condition_groups[[i]],
                        pbm_conditions_group = names(pbm_condition_groups[i]),
                        check_corecmotifs = FALSE
                    )
                }) %>%

                # Combine all the data frames
                dplyr::bind_rows()

            try_catch_save_output(
                condition_comparison_df, comparison_output, "tsv"
            )
        }
    }

    # Return the list of CoRecMotifs
    return(final_corecmotifs)
}

