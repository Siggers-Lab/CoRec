#' Filter a list of \linkS4class{CoRecMotif} objects
#'
#' Filter a list of \linkS4class{CoRecMotif} objects based on the contents of
#' \linkS4class{CoRecMotif} slots.
#'
#' @param corecmotifs the list of CoRecMotif objects to filter.
#' @param seed_name a character vector specifying the seed names (probe sets) to
#'   keep.
#' @param pbm_condition a character vector specifying the PBM conditions to
#'   keep.
#' @param array_id a character vector specifying the array IDs to keep.
#' @param rolling_ic a single number specifying the minimum rolling IC to keep.
#' @param motif_strength a single number specifying the minimum motif strength
#'   to keep.
#' @param seed_sequence a character vector specifying the seed sequences to
#'   keep.
#' @param motif_name a character vector specifying the motif names to keep.
#' @param match_name a character vector specifying the matching motif names to
#'   keep.
#' @param match_altname a character vector specifying the matching motif
#'   altnames to keep.
#' @param match_pvalue a single number specifying the maximum match p-value to
#'   keep.
#' @param match_cluster a character vector specifying the matching clusters to
#'   keep.
#'
#' @return A list of \linkS4class{CoRecMotif} objects that pass the filters.
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
#' # Filter by rolling IC and motif strength
#' filtered_corecmotifs <-
#'     filter_corecmotifs(
#'         corecmotifs,
#'         rolling_ic = 1,
#'         motif_strength = 1
#'     )
filter_corecmotifs <-
    function(
        corecmotifs,
        seed_name = NULL,
        pbm_condition = NULL,
        array_id = NULL,
        rolling_ic = NULL,
        motif_strength = NULL,
        seed_sequence = NULL,
        motif_name = NULL,
        match_name = NULL,
        match_altname = NULL,
        match_pvalue = NULL,
        match_cluster = NULL
    ) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        is.character(seed_name) || is.null(seed_name),
        msg = "seed_name is not a character vector or NULL"
    )
    assertthat::assert_that(
        is.character(pbm_condition) || is.null(pbm_condition),
        msg = "pbm_condition is not a character vector or NULL"
    )
    assertthat::assert_that(
        is.character(array_id) || is.null(array_id),
        msg = "array_id is not a character vector or NULL"
    )
    assertthat::assert_that(
        assertthat::is.number(rolling_ic) || is.null(rolling_ic),
        msg = "rolling_ic is not a single number or NULL"
    )
    assertthat::assert_that(
        assertthat::is.number(motif_strength) || is.null(motif_strength),
        msg = "motif_strength is not a single number or NULL"
    )
    assertthat::assert_that(
        is.character(seed_sequence) || is.null(seed_sequence),
        msg = "seed_sequence is not a character vector or NULL"
    )
    assertthat::assert_that(
        is.character(motif_name) || is.null(motif_name),
        msg = "motif_name is not a character vector or NULL"
    )
    assertthat::assert_that(
        is.character(match_name) || is.null(match_name),
        msg = "match_name is not a character vector or NULL"
    )
    assertthat::assert_that(
        is.character(match_altname) || is.null(match_altname),
        msg = "match_altname is not a character vector or NULL"
    )
    assertthat::assert_that(
        assertthat::is.number(match_pvalue) || is.null(match_pvalue),
        msg = "match_pvalue is not a single number or NULL"
    )
    assertthat::assert_that(
        is.character(match_cluster) || is.null(match_cluster),
        msg = "match_cluster is not a character vector or NULL"
    )

    # Make sure corecmotifs is a list
    if (!is.list(corecmotifs)) {
        corecmotifs <- list(corecmotifs)
    }

    # Filter by seed_name
    if (!is.null(seed_name)) {
        vals <- vapply(corecmotifs, get_seed_name, character(1))
        corecmotifs <- corecmotifs[vals %in% seed_name]
    }

    # Filter by pbm_condition
    if (!is.null(pbm_condition)) {
        vals <- vapply(corecmotifs, get_pbm_condition, character(1))
        corecmotifs <- corecmotifs[vals %in% pbm_condition]
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

    # Return the filtered list
    return(corecmotifs)
}

