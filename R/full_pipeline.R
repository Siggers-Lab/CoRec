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
        output_directory = NULL,
        output_base_name = NULL,
        array_id = NULL
    ) {
    # Make sure the output arguments are the right type
    # Everything else is checked in the other functions this one calls
    assertthat::assert_that(
        assertthat::is.string(output_directory) || is.null(output_directory),
        msg = "output_directory is not a character vector or NULL"
    )
    assertthat::assert_that(
        assertthat::is.string(output_base_name) || is.null(output_base_name),
        msg = "output_base_name is not a character vector or NULL"
    )
    assertthat::assert_that(
        assertthat::is.string(array_id) || is.null(array_id),
        msg = "array_id is not a character vector or NULL"
    )

    # Do not save any output files by default
    fluorescence_output <- NULL
    zscore_output <- NULL
    corec_output <- NULL

    # If only the output directory or base name is provided, set the other
    if (is.null(output_directory) & !is.null(output_base_name)) {
        # If no output directory is provided, use current working directory
        output_directory <- getwd()
    } else if (!is.null(output_directory) & is.null(output_base_name)) {
        # If no output base name is provided, use "output"
        output_base_name <- "output"
    }

    # If both the output directory and base name exist, set output file names
    if (!is.null(output_directory) & !is.null(output_base_name)) {
        # Create the output directory if it doesn't already exist
        if (!dir.exists(output_directory)) {
            dir.create(output_directory)
        }

        # Add the output directory path to the base file name for output files
        output_base_name <- paste(output_directory, output_base_name, sep = "/")

        # Add the array ID to the base file name if it's provided
        if (!is.null(array_id)) {
            output_base_name <- paste(output_base_name, array_id, sep = "_")
        }

        # Create names for the three output files
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

    # Make corecmotif objects for all the seed/condition combos
    corecmotifs <-
        zscore_table_to_corecmotifs(
            zscore_table = zscore_table,
            zscore_columns = pbm_conditions,
            output_file = corec_output,
            array_id = array_id
        )

    # Return the list of corecmotifs
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
    ) {
    # Filter out corecmotifs with low motif strength and/or rolling IC scores
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
            meme_path = meme_path
        )

    # Filter out CoRecMotifs that don't match any reference motifs well
    filtered_matched_corecmotifs <-
        filter_corecmotifs(
            matched_corecmotifs,
            match_pvalue = max_match_pvalue
        ) %>%

        # Make sure at least min_n_replicates match a reference motif well
        match_replicates(
            min_n_replicates = min_n_replicates,
            max_eucl_distance = NULL
        )
}

