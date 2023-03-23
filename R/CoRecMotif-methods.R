#' @param x,object [CoRecMotif][CoRecMotif-class]. The motif to manipulate.
#'
#' @export
#' @rdname CoRecMotif-class
methods::setMethod("as.data.frame", "CoRecMotif", function(x) {
    data.frame(
        probe_set = get_probe_set(x),
        pbm_condition = get_pbm_condition(x),
        array_id = get_array_id(x),
        motif_strength = get_motif_strength(x),
        rolling_ic = get_rolling_ic(x),
        seed_sequence = get_seed_sequence(x),
        match_motif = ifelse(
            methods::is(get_match_motif(x), "universalmotif"),
            get_match_altname(x),
            NA
        ),
        match_pvalue = get_match_pvalue(x),
        match_cluster = get_match_cluster(x),
        check.names = FALSE,
        fix.empty.names = FALSE,
        stringsAsFactors = FALSE
    )
})

#' @export
#' @rdname CoRecMotif-class
methods::setMethod("show", "CoRecMotif", function(object) {
    cat("\n",
        "   Probe set:       ", get_probe_set(object), "\n",
        "   PBM condition:   ", get_pbm_condition(object), "\n",
        "   Array ID:        ", get_array_id(object), "\n",
        "   Motif strength:  ", round(get_motif_strength(object), 2), "\n",
        "   Rolling IC:      ", round(get_rolling_ic(object), 2), "\n",
        "   Match motif:     ", get_match_altname(object), "\n",
        "   Match cluster:   ", get_match_cluster(object), "\n\n",
        "Delta z-score motif:\n\n",
        sep = ""
    )
    print(round(get_delta_zscore_motif(object), 1))
    cat("\nPPM:\n\n")
    print(round(get_ppm(object), 2))
})

#' Access or replace data from CoRecMotifs
#'
#' Extracts or replaces data stored in [CoRecMotifs][CoRecMotif-class].
#'
#' @param corecmotif [CoRecMotif][CoRecMotif-class]. The motif whose slots
#'   should be accessed or replaced.
#' @param correct_orientation `logical(1)`. Should the reference motif be
#'   reversed if necessary to match the CoRecMotif's orientation? (Default:
#'   TRUE)
#' @param value Object to replace the slot with. The required type depends on
#'   which slot is being replaced. See [CoRecMotif-class] for more information
#'   on the expected classes.
#'
#' @return The value of the selected slot OR a [CoRecMotif][CoRecMotif-class]
#'   with the appropriate slot updated.
#'
#' @export
#' @name accessors
#' @rdname accessors
get_probe_set <- function(corecmotif) {
    return(corecmotif@probe_set)
}

#' @export
#' @rdname accessors
get_pbm_condition <- function(corecmotif) {
    return(corecmotif@pbm_condition)
}

#' @export
#' @rdname accessors
get_zscore_motif <- function(corecmotif) {
    return(corecmotif@zscore_motif)
}

#' @export
#' @rdname accessors
get_delta_zscore_motif <- function(corecmotif) {
    # Get the z-score motif
    zscore_motif <- get_zscore_motif(corecmotif)

    # Transform the z-scores to reflect their deviation from column-wise median
    delta_zscore_motif <-
        # Subtract the column-wise median from each value in each column
        apply(zscore_motif, 2, function(col) col - stats::median(col))

    # Return the delta z-score motif
    return(delta_zscore_motif)
}

#' @export
#' @rdname accessors
get_array_id <- function(corecmotif) {
    return(corecmotif@array_id)
}

#' @export
#' @rdname accessors
get_motif_strength <- function(corecmotif) {
    return(corecmotif@motif_strength)
}

#' @export
#' @rdname accessors
get_rolling_ic <- function(corecmotif) {
    return(corecmotif@rolling_ic)
}

#' @export
#' @rdname accessors
get_seed_zscore <- function(corecmotif) {
    # Get the z-score motif
    zscore_motif <- get_zscore_motif(corecmotif)

    # Find the seed probe z-score
    seed_zscore <- find_seed_zscore(zscore_motif)

    # Return the z-score
    return(seed_zscore)
}

#' @export
#' @rdname accessors
get_seed_sequence <- function(corecmotif) {
    return(corecmotif@seed_sequence)
}

#' @export
#' @rdname accessors
get_beta <- function(corecmotif) {
    return(corecmotif@beta)
}

#' @export
#' @rdname accessors
get_motif <- function(corecmotif) {
    return(corecmotif@motif)
}

#' @export
#' @rdname accessors
get_motif_name <- function(corecmotif) {
    return(corecmotif@motif@name)
}

#' @export
#' @rdname accessors
get_icm <- function(corecmotif) {
    icm <- universalmotif::convert_type(corecmotif@motif, "ICM")
    return(icm["motif"])
}

#' @export
#' @rdname accessors
get_pwm <- function(corecmotif) {
    pwm <- universalmotif::convert_type(corecmotif@motif, "PWM")
    return(pwm["motif"])
}

#' @export
#' @rdname accessors
get_ppm <- function(corecmotif) {
    ppm <- universalmotif::convert_type(corecmotif@motif, "PPM")
    return(ppm["motif"])
}

#' @export
#' @rdname accessors
get_match_motif <- function(corecmotif) {
    return(corecmotif@match_motif)
}

#' @export
#' @rdname accessors
get_match_name <- function(corecmotif) {
    if (methods::is(corecmotif@match_motif, "universalmotif")) {
        return(corecmotif@match_motif@name)
    } else {
        return(NA)
    }
}

#' @export
#' @rdname accessors
get_match_altname <- function(corecmotif) {
    if (methods::is(corecmotif@match_motif, "universalmotif")) {
        return(corecmotif@match_motif@altname)
    } else {
        return(NA)
    }
}

#' @export
#' @rdname accessors
get_match_icm <- function(corecmotif, correct_orientation = TRUE) {
    if (methods::is(corecmotif@match_motif, "universalmotif")) {
        icm <- universalmotif::convert_type(corecmotif@match_motif, "ICM")
        if (correct_orientation && icm["extrainfo"] == "-") {
            icm <- universalmotif::motif_rc(icm)
        }
        return(icm["motif"])
    } else {
        return(NA)
    }
}

#' @export
#' @rdname accessors
get_match_pwm <- function(corecmotif, correct_orientation = TRUE) {
    if (methods::is(corecmotif@match_motif, "universalmotif")) {
        pwm <- universalmotif::convert_type(corecmotif@match_motif, "PWM")
        if (correct_orientation && pwm["extrainfo"] == "-") {
            pwm <- universalmotif::motif_rc(pwm)
        }
        return(pwm["motif"])
    } else {
        return(NA)
    }
}

#' @export
#' @rdname accessors
get_match_ppm <- function(corecmotif, correct_orientation = TRUE) {
    if (methods::is(corecmotif@match_motif, "universalmotif")) {
        ppm <- universalmotif::convert_type(corecmotif@match_motif, "PPM")
        if (correct_orientation && ppm["extrainfo"] == "-") {
            ppm <- universalmotif::motif_rc(ppm)
        }
        return(ppm["motif"])
    } else {
        return(NA)
    }
}

#' @export
#' @rdname accessors
get_match_pvalue <- function(corecmotif) {
    return(corecmotif@match_pvalue)
}

#' @export
#' @rdname accessors
get_match_cluster <- function(corecmotif) {
    return(corecmotif@match_cluster)
}

#' @export
#' @rdname accessors
set_probe_set <- function(corecmotif, value) {
    corecmotif@probe_set <- value
    methods::validObject(corecmotif)
    return(corecmotif)
}

#' @export
#' @rdname accessors
set_pbm_condition <- function(corecmotif, value) {
    corecmotif@pbm_condition <- value
    methods::validObject(corecmotif)
    return(corecmotif)
}

#' @export
#' @rdname accessors
set_zscore_motif <- function(corecmotif, value) {
    # Make sure zscore_motif is the right format
    value <- check_valid_zscore_motif(value)

    # Update the z-score motif itself
    corecmotif@zscore_motif <- value

    # Update all the things that depend on the z-score motif
    corecmotif@motif_strength <- calculate_strength(value)
    corecmotif@beta <- calculate_beta(corecmotif@motif_strength)
    corecmotif@motif <- zscore_to_universalmotif(
        value, corecmotif@beta, get_motif_name(corecmotif)
    )
    corecmotif@rolling_ic <- calculate_rolling_ic(corecmotif@motif)

    # Reset all the match slots
    corecmotif@match_motif <- NA
    corecmotif@match_pvalue <- NA_real_
    corecmotif@match_cluster <- NA_character_

    methods::validObject(corecmotif)
    return(corecmotif)
}

#' @export
#' @rdname accessors
set_array_id <- function(corecmotif, value) {
    corecmotif@array_id <- value
    methods::validObject(corecmotif)
    return(corecmotif)
}

#' @export
#' @rdname accessors
set_seed_sequence <- function(corecmotif, value) {
    corecmotif@seed_sequence <- value
    methods::validObject(corecmotif)
    return(corecmotif)
}

#' @export
#' @rdname accessors
set_motif_name <- function(corecmotif, value) {
    corecmotif@motif@name <- value
    methods::validObject(corecmotif)
    return(corecmotif)
}

