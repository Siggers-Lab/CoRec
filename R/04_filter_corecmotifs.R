#' Filter a list of CoRecMotifs
#'
#' Filter a list of [CoRecMotifs][CoRecMotif-class] based on the contents of
#' CoRecMotif slots.
#'
#' To skip filtering by a particular slot, set the corresponding argument to
#' NULL (the default).
#'
#' @param corecmotifs the list of [CoRecMotifs][CoRecMotif-class] to filter.
#' @param probe_set a character vector specifying the probe sets to keep.
#'   (Default: NULL)
#' @param pbm_condition a character vector specifying the PBM conditions to
#'   keep. (Default: NULL)
#' @param array_id a character vector specifying the array IDs to keep.
#'   (Default: NULL)
#' @param rolling_ic a single number specifying the minimum rolling IC to keep.
#'   (Default: NULL)
#' @param motif_strength a single number specifying the minimum motif strength
#'   to keep. (Default: NULL)
#' @param seed_sequence a character vector specifying the seed sequences to
#'   keep. (Default: NULL)
#' @param motif_name a character vector specifying the motif names to keep.
#'   (Default: NULL)
#' @param match_name a character vector specifying the matching motif names to
#'   keep. (Default: NULL)
#' @param match_altname a character vector specifying the matching motif
#'   altnames to keep. (Default: NULL)
#' @param match_pvalue a single number specifying the maximum match p-value to
#'   keep. (Default: NULL)
#' @param match_cluster a character vector specifying the matching clusters to
#'   keep. (Default: NULL)
#' @param output_file the path to the RDS file where the filtered list of
#'   [CoRecMotifs][CoRecMotif-class] will be written. If NULL, no file is
#'   written. (Default: NULL)
#'
#' @return A list of [CoRecMotifs][CoRecMotif-class] that pass the filters.
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN")
filter_corecmotifs <-
    function(
        corecmotifs,
        probe_set = NULL,
        pbm_condition = NULL,
        array_id = NULL,
        rolling_ic = NULL,
        motif_strength = NULL,
        seed_sequence = NULL,
        motif_name = NULL,
        match_name = NULL,
        match_altname = NULL,
        match_pvalue = NULL,
        match_cluster = NULL,
        output_file = NULL
    ) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        is.character(probe_set) || is.null(probe_set),
        is.character(pbm_condition) || is.null(pbm_condition),
        is.character(array_id) || is.null(array_id),
        assertthat::is.number(rolling_ic) || is.null(rolling_ic),
        assertthat::is.number(motif_strength) || is.null(motif_strength),
        is.character(seed_sequence) || is.null(seed_sequence),
        is.character(motif_name) || is.null(motif_name),
        is.character(match_name) || is.null(match_name),
        is.character(match_altname) || is.null(match_altname),
        assertthat::is.number(match_pvalue) || is.null(match_pvalue),
        is.character(match_cluster) || is.null(match_cluster),
        assertthat::is.string(output_file) || is.null(output_file)
    )

    # Make sure corecmotifs is a list
    if (!is.list(corecmotifs)) {
        corecmotifs <- list(corecmotifs)
    }

    # Filter by probe_set
    if (!is.null(probe_set)) {
        vals <- vapply(corecmotifs, get_probe_set, character(1))
        corecmotifs <- corecmotifs[vals %in% probe_set]
    }

    # Filter by pbm_condition
    if (!is.null(pbm_condition)) {
        vals <- vapply(corecmotifs, get_pbm_condition, character(1))
        corecmotifs <- corecmotifs[vals %in% pbm_condition]
    }

    # Filter by array_id
    if (!is.null(array_id)) {
        vals <- vapply(corecmotifs, get_array_id, character(1))
        corecmotifs <- corecmotifs[vals %in% array_id]
    }

    # Filter by rolling_ic
    if (!is.null(rolling_ic)) {
        vals <- vapply(corecmotifs, get_rolling_ic, numeric(1))
        corecmotifs <- corecmotifs[vals >= rolling_ic & !is.na(vals)]
    }

    # Filter by motif_strength
    if (!is.null(motif_strength)) {
        vals <- vapply(corecmotifs, get_motif_strength, numeric(1))
        corecmotifs <- corecmotifs[vals >= motif_strength & !is.na(vals)]
    }

    # Filter by seed_sequence
    if (!is.null(seed_sequence)) {
        vals <- vapply(corecmotifs, get_seed_sequence, character(1))
        corecmotifs <- corecmotifs[vals %in% seed_sequence]
    }

    # Filter by motif_name
    if (!is.null(motif_name)) {
        vals <- vapply(corecmotifs, get_motif_name, character(1))
        corecmotifs <- corecmotifs[vals %in% motif_name]
    }

    # Filter by match_name
    if (!is.null(match_name)) {
        vals <- vapply(corecmotifs, get_match_name, character(1))
        corecmotifs <- corecmotifs[vals %in% match_name]
    }

    # Filter by match_altname
    if (!is.null(match_altname)) {
        vals <- vapply(corecmotifs, get_match_altname, character(1))
        corecmotifs <- corecmotifs[vals %in% match_altname]
    }

    # Filter by match_pvalue
    if (!is.null(match_pvalue)) {
        vals <- vapply(corecmotifs, get_match_pvalue, numeric(1))
        corecmotifs <- corecmotifs[vals <= match_pvalue & !is.na(vals)]
    }

    # Filter by match_cluster
    if (!is.null(match_cluster)) {
        vals <- vapply(corecmotifs, get_match_cluster, character(1))
        corecmotifs <- corecmotifs[vals %in% match_cluster]
    }

    # Try to save the filtered CoRecMotifs as an RDS file if necessary
    try_catch_save_output(corecmotifs, output_file, "rds")

    # Return the filtered list
    return(corecmotifs)
}

