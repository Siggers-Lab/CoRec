# Create generics for all the necessary getters --------------------------------
setGeneric("probe_set", function(x) standardGeneric("probe_set"))
setGeneric("pbm_condition", function(x) standardGeneric("pbm_condition"))
setGeneric("zscore_motif", function(x) standardGeneric("zscore_motif"))
setGeneric(
    "delta_zscore_motif", function(x) standardGeneric("delta_zscore_motif")
)
setGeneric("array_id", function(x) standardGeneric("array_id"))
setGeneric("motif_strength", function(x) standardGeneric("motif_strength"))
setGeneric("rolling_ic", function(x) standardGeneric("rolling_ic"))
setGeneric("seed_zscore", function(x) standardGeneric("seed_zscore"))
setGeneric("seed_sequence", function(x) standardGeneric("seed_sequence"))
setGeneric("beta", function(x) standardGeneric("beta"))

setGeneric("motif", function(x) standardGeneric("motif"))
setGeneric("motif_name", function(x) standardGeneric("motif_name"))
setGeneric("icm", function(x) standardGeneric("icm"))
setGeneric("pwm", function(x) standardGeneric("pwm"))
setGeneric("ppm", function(x) standardGeneric("ppm"))

setGeneric("match_motif", function(x) standardGeneric("match_motif"))
setGeneric("match_name", function(x) standardGeneric("match_name"))
setGeneric("match_altname", function(x) standardGeneric("match_altname"))
setGeneric("match_icm", function(x) standardGeneric("match_icm"))
setGeneric("match_pwm", function(x) standardGeneric("match_pwm"))
setGeneric("match_ppm", function(x) standardGeneric("match_ppm"))
setGeneric("match_pvalue", function(x) standardGeneric("match_pvalue"))
setGeneric("match_cluster", function(x) standardGeneric("match_cluster"))

# Define the methods for all the getters

#' @rdname CoRecMotif-class
setMethod("probe_set", "CoRecMotif", function(x) x@probe_set)

#' @rdname CoRecMotif-class
setMethod("pbm_condition", "CoRecMotif", function(x) x@pbm_condition)

#' @rdname CoRecMotif-class
setMethod("zscore_motif", "CoRecMotif", function(x) x@zscore_motif)

#' @rdname CoRecMotif-class
setMethod("delta_zscore_motif", "CoRecMotif", function(x) {
    # Get the z-score motif
    zscore_motif <- zscore_motif(x)

    # Transform the z-scores to reflect their deviation from column-wise median
    delta_zscore_motif <-
        # Subtract the column-wise median from each value in each column
        apply(zscore_motif, 2, function(col) col - median(col))

    # Return the delta z-score motif
    return(delta_zscore_motif)
})

#' @rdname CoRecMotif-class
setMethod("array_id", "CoRecMotif", function(x) x@array_id)

#' @rdname CoRecMotif-class
setMethod("motif_strength", "CoRecMotif", function(x) x@motif_strength)

#' @rdname CoRecMotif-class
setMethod("rolling_ic", "CoRecMotif", function(x) x@rolling_ic)

#' @rdname CoRecMotif-class
setMethod("seed_zscore", "CoRecMotif", function(x) {
    # Get the z-score motif
    zscore_motif <- zscore_motif(x)

    # Find the seed probe z-score
    seed_zscore <- find_seed_zscore(zscore_motif)

    # Return the z-score
    return(seed_zscore)
})

#' @rdname CoRecMotif-class
setMethod("seed_sequence", "CoRecMotif", function(x) x@seed_sequence)

#' @rdname CoRecMotif-class
setMethod("beta", "CoRecMotif", function(x) x@beta)

# @motif related getters

#' @rdname CoRecMotif-class
setMethod("motif", "CoRecMotif", function(x) x@motif)

#' @rdname CoRecMotif-class
setMethod("motif_name", "CoRecMotif", function(x) x@motif@name)

#' @rdname CoRecMotif-class
setMethod("icm", "CoRecMotif", function(x) {
    universalmotif::convert_type(x@motif, "ICM")@motif
})

#' @rdname CoRecMotif-class
setMethod("pwm", "CoRecMotif", function(x) {
    universalmotif::convert_type(x@motif, "PWM")@motif
})

#' @rdname CoRecMotif-class
setMethod("ppm", "CoRecMotif", function(x) {
    universalmotif::convert_type(x@motif, "PPM")@motif
})

# @match_motif related getters

#' @rdname CoRecMotif-class
setMethod("match_motif", "CoRecMotif", function(x) x@match_motif)

#' @rdname CoRecMotif-class
setMethod("match_name", "CoRecMotif", function(x) {
    if (is(x@match_motif, "universalmotif")) {
        x@match_motif@name
    } else {
        NA
    }
})

#' @rdname CoRecMotif-class
setMethod("match_altname", "CoRecMotif", function(x) {
    if (is(x@match_motif, "universalmotif")) {
        x@match_motif@altname
    } else {
        NA
    }
})

#' @rdname CoRecMotif-class
setMethod("match_icm", "CoRecMotif", function(x) {
    if (is(x@match_motif, "universalmotif")) {
        universalmotif::convert_type(x@match_motif, "ICM")@motif
    } else {
        NA
    }
})

#' @rdname CoRecMotif-class
setMethod("match_pwm", "CoRecMotif", function(x) {
    if (is(x@match_motif, "universalmotif")) {
        universalmotif::convert_type(x@match_motif, "PWM")@motif
    } else {
        NA
    }
})

#' @rdname CoRecMotif-class
setMethod("match_ppm", "CoRecMotif", function(x) {
    if (is(x@match_motif, "universalmotif")) {
        universalmotif::convert_type(x@match_motif, "PPM")@motif
    } else {
        NA
    }
})

#' @rdname CoRecMotif-class
setMethod("match_pvalue", "CoRecMotif", function(x) x@match_pvalue)

#' @rdname CoRecMotif-class
setMethod("match_cluster", "CoRecMotif", function(x) x@match_cluster)

# Create generics for all the necessary setters --------------------------------
setGeneric(
    "probe_set<-", function(x, value) standardGeneric("probe_set<-")
)
setGeneric(
    "pbm_condition<-", function(x, value) standardGeneric("pbm_condition<-")
)
setGeneric(
    "zscore_motif<-", function(x, value) standardGeneric("zscore_motif<-")
)
setGeneric(
    "array_id<-", function(x, value) standardGeneric("array_id<-")
)
setGeneric(
    "seed_sequence<-", function(x, value) standardGeneric("seed_sequence<-")
)
setGeneric(
    "motif_name<-", function(x, value) standardGeneric("motif_name<-")
)

# Define the methods for all the setters

#' @rdname CoRecMotif-class
setMethod("probe_set<-", "CoRecMotif", function(x, value) {
    x@probe_set <- value
    validObject(x)
    x
})

#' @rdname CoRecMotif-class
setMethod("pbm_condition<-", "CoRecMotif", function(x, value) {
    x@pbm_condition <- value
    validObject(x)
    x
})

#' @rdname CoRecMotif-class
setMethod("zscore_motif<-", "CoRecMotif", function(x, value) {
    # Check the input to start with or you'll get unhelpful error messages later
    if (!is.numeric(value) || !is.matrix(value)) {
        stop(
            "zscore_motif must be a numeric matrix",
            call. = FALSE
        )
    }

    # Update the z-score motif itself
    x@zscore_motif <- value

    # Update all the things that depend on the z-score motif
    x@motif_strength <- calculate_strength(value)
    x@beta <- calculate_beta(value)
    x@motif <- zscore_to_universalmotif(value, x@beta, motif_name(x))
    x@rolling_ic <- calculate_rolling_ic(x@motif)

    # Reset all the match slots
    x@match_motif <- NA
    x@match_pvalue <- NA_real_
    x@match_cluster <- NA_character_

    validObject(x)
    x
})

#' @rdname CoRecMotif-class
setMethod("array_id<-", "CoRecMotif", function(x, value) {
    x@array_id <- value
    validObject(x)
    x
})

#' @rdname CoRecMotif-class
setMethod("seed_sequence<-", "CoRecMotif", function(x, value) {
    x@seed_sequence <- value
    validObject(x)
    x
})

#' @rdname CoRecMotif-class
setMethod("motif_name<-", "CoRecMotif", function(x, value) {
    x@motif@name <- value
    validObject(x)
    x
})

# Define the methods for other useful generics ---------------------------------

#' @rdname CoRecMotif-class
setMethod("as.data.frame", "CoRecMotif", function(x) {
    data.frame(
        probe_set = probe_set(x),
        pbm_condition = pbm_condition(x),
        array_id = array_id(x),
        rolling_ic = rolling_ic(x),
        motif_strength = motif_strength(x),
        seed_sequence = seed_sequence(x),
        match_motif = ifelse(
            is(match_motif(x), "universalmotif"),
            match_altname(x),
            NA
        ),
        match_pvalue = match_pvalue(x),
        match_cluster = match_cluster(x),
        check.names = FALSE,
        fix.empty.names = FALSE,
        stringsAsFactors = FALSE
    )
})

#' @rdname CoRecMotif-class
setMethod("show", "CoRecMotif", function(object) {
    cat("\n",
        "   Probe set:      ", probe_set(object), "\n",
        "   PBM condition:  ", pbm_condition(object), "\n",
        "   Array ID:       ", array_id(object), "\n",
        "   Motif strength: ", round(motif_strength(object), 2), "\n",
        "   Rolling IC:     ", round(rolling_ic(object), 2), "\n\n"
    )
    print(round(ppm(object), 2))
})

