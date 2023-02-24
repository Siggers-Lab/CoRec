#' Filter CoRecMotifs
#'
#' Filter a list of [CoRecMotifs][CoRecMotif-class], retaining only those that
#' satisfy all of your conditions.
#'
#' @inheritParams annotate_fluorescence_table
#' @param corecmotifs `list`. The [CoRecMotifs][CoRecMotif-class] to filter.
#' @param probe_set `character` or `NULL`. The probe sets to keep or NULL not to
#'   filter by probe set. (Default: NULL)
#' @param pbm_condition `character` or `NULL`. The PBM conditions to keep or
#'   NULL not to filter by PBM condition. (Default: NULL)
#' @param array_id `character` or `NULL`. The array IDs to keep or NULL not to
#'   filter by array ID. (Default: NULL)
#' @param rolling_ic `numeric(1)` or `NULL`. The minimum rolling IC to keep or
#'   NULL not to filter by rolling IC. (Default: NULL)
#' @param motif_strength `numeric(1)` or `NULL`. The minimum motif strength to
#'   keep or NULL not to filter by motif strength. (Default: NULL)
#' @param seed_sequence `character` or `NULL`. The seed sequences to keep or
#'   NULL not to filter by seed sequence. (Default: NULL)
#' @param motif_name `character` or `NULL`. The motif names to keep or NULL not
#'   to filter by motif name. (Default: NULL)
#' @param match_name `character` or `NULL`. The match motif names to keep or
#'   NULL not to filter by match motif name. (Default: NULL)
#' @param match_altname `character` or `NULL`. The match motif altnames to keep
#'   or NULL not to filter by match motif altname. (Default: NULL)
#' @param match_pvalue `numeric(1)` or `NULL`. The maximum match p-value to keep
#'   or NULL not to filter by match p-value. (Default: NULL)
#' @param match_cluster `character` or `NULL`. The match clusters to keep or
#'   NULL not to filter by match cluster. (Default: NULL)
#'
#' @return A list of [CoRecMotifs][CoRecMotif-class] that pass the given
#'   filters.
#'
#' @seealso [check_replicates()] for filtering groups of replicate
#'   [CoRecMotifs][CoRecMotif-class].
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
        motif_strength = NULL,
        rolling_ic = NULL,
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
        assertthat::is.number(motif_strength) || is.null(motif_strength),
        assertthat::is.number(rolling_ic) || is.null(rolling_ic),
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

    # Filter by motif_strength
    if (!is.null(motif_strength)) {
        vals <- vapply(corecmotifs, get_motif_strength, numeric(1))
        corecmotifs <- corecmotifs[vals >= motif_strength & !is.na(vals)]
    }

    # Filter by rolling_ic
    if (!is.null(rolling_ic)) {
        vals <- vapply(corecmotifs, get_rolling_ic, numeric(1))
        corecmotifs <- corecmotifs[vals >= rolling_ic & !is.na(vals)]
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

