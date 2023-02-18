setMethod(
    "as.data.frame",
    signature(x = "corecmotif"),
    definition = function(x) {
        data.frame(
            seed_name = get_seed_name(x),
            pbm_condition = get_pbm_condition(x),
            array_id = get_array_id(x),
            rolling_ic = get_rolling_ic(x),
            motif_strength = get_motif_strength(x),
            seed_sequence = get_seed_sequence(x),
            match_motif = ifelse(
                is.null(x@match_motif),
                NA,
                x@match_motif@altname
            ),
            match_pvalue = get_match_pvalue(x),
            match_cluster = get_match_cluster(x),
            check.names = FALSE,
            fix.empty.names = FALSE,
            stringsAsFactors = FALSE
        )
    })
