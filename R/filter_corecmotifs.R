filter_corecmotifs <-
    function(
        corecmotifs,
        seed_name = NULL,
        pbm_condition = NULL,
        rolling_ic = NULL,
        motif_strength = NULL,
        seed_sequence = NULL,
        match_method = NULL,
        match_pvalue = NULL,
        match_cluster = NULL
    ) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        is.character(seed_name) | is.null(seed_name),
        msg = "seed_name must be a character vector or NULL"
    )
    assertthat::assert_that(
        is.character(pbm_condition) | is.null(pbm_condition),
        msg = "pbm_condition must be a character vector or NULL"
    )
    assertthat::assert_that(
        assertthat::is.number(rolling_ic) | is.null(rolling_ic),
        msg = "rolling_ic must be a single number or NULL"
    )
    assertthat::assert_that(
        assertthat::is.number(motif_strength) | is.null(motif_strength),
        msg = "motif_strength must be a single number or NULL"
    )
    assertthat::assert_that(
        is.character(seed_sequence) | is.null(seed_sequence),
        msg = "seed_sequence must be a character vector or NULL"
    )
    assertthat::assert_that(
        is.character(match_method) | is.null(match_method),
        msg = "match_method must be a character vector or NULL"
    )
    assertthat::assert_that(
        assertthat::is.number(match_pvalue) | is.null(match_pvalue),
        msg = "match_pvalue must be a single number or NULL"
    )
    assertthat::assert_that(
        is.character(match_cluster) | is.null(match_cluster),
        msg = "match_cluster must be a character vector or NULL"
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

    # Filter by match_method
    if (!is.null(match_method)) {
        vals <- vapply(corecmotifs, get_match_method, character(1))
        corecmotifs <- corecmotifs[vals %in% match_method]
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

