# Create generics for all the necessary getters --------------------------------
setGeneric("probe_set", function(x) standardGeneric("probe_set"))
setGeneric("pbm_condition", function(x) standardGeneric("pbm_condition"))
setGeneric("zscore_motif", function(x) standardGeneric("zscore_motif"))
setGeneric("array_id", function(x) standardGeneric("array_id"))
setGeneric("motif_strength", function(x) standardGeneric("motif_strength"))
setGeneric("rolling_ic", function(x) standardGeneric("rolling_ic"))
setGeneric("seed_sequence", function(x) standardGeneric("seed_sequence"))
setGeneric("beta", function(x) standardGeneric("beta"))
setGeneric("ppm", function(x) standardGeneric("ppm"))
setGeneric("motif_name", function(x) standardGeneric("motif_name"))
setGeneric("match_motif", function(x) standardGeneric("match_motif"))
setGeneric("match_name", function(x) standardGeneric("match_name"))
setGeneric("match_altname", function(x) standardGeneric("match_altname"))
setGeneric("match_pvalue", function(x) standardGeneric("match_pvalue"))
setGeneric("match_cluster", function(x) standardGeneric("match_cluster"))
setGeneric("seed_zscore", function(x) standardGeneric("seed_zscore"))
setGeneric(
    "delta_zscore_motif", function(x) standardGeneric("delta_zscore_motif")
)

# Define the methods for all the getters

#' @rdname CoRecMotif-class
setMethod("probe_set", "CoRecMotif", function(x) x@probe_set)

#' @rdname CoRecMotif-class
setMethod("pbm_condition", "CoRecMotif", function(x) x@pbm_condition)

#' @rdname CoRecMotif-class
setMethod("zscore_motif", "CoRecMotif", function(x) x@zscore_motif)

#' @rdname CoRecMotif-class
setMethod("array_id", "CoRecMotif", function(x) x@array_id)

#' @rdname CoRecMotif-class
setMethod("motif_strength", "CoRecMotif", function(x) x@motif_strength)

#' @rdname CoRecMotif-class
setMethod("rolling_ic", "CoRecMotif", function(x) x@rolling_ic)

#' @rdname CoRecMotif-class
setMethod("seed_sequence", "CoRecMotif", function(x) x@seed_sequence)

#' @rdname CoRecMotif-class
setMethod("beta", "CoRecMotif", function(x) x@beta)

#' @rdname CoRecMotif-class
setMethod("ppm", "CoRecMotif", function(x) x@ppm)

#' @rdname CoRecMotif-class
setMethod("motif_name", "CoRecMotif", function(x) x@ppm@name)

#' @rdname CoRecMotif-class
setMethod("match_motif", "CoRecMotif", function(x) x@match_motif)

#' @rdname CoRecMotif-class
setMethod("match_name", "CoRecMotif", function(x) {
    if (is(x@ppm, "universalmotif")) {
        x@match_motif@name
    } else {
        NA
    }
})

#' @rdname CoRecMotif-class
setMethod("match_altname", "CoRecMotif", function(x) {
    if (is(x@ppm, "universalmotif")) {
        x@match_motif@altname
    } else {
        NA
    }
})

#' @rdname CoRecMotif-class
setMethod("match_pvalue", "CoRecMotif", function(x) x@match_pvalue)

#' @rdname CoRecMotif-class
setMethod("match_cluster", "CoRecMotif", function(x) x@match_cluster)

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
setGeneric(
    "match_motif<-", function(x, value) standardGeneric("match_motif<-")
)
setGeneric(
    "match_pvalue<-", function(x, value) standardGeneric("match_pvalue<-")
)
setGeneric(
    "match_cluster<-", function(x, value) standardGeneric("match_cluster<-")
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
    x@zscore_motif <- value
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
    x@ppm@name <- value
    validObject(x)
    x
})

#' @rdname CoRecMotif-class
setMethod("match_motif<-", "CoRecMotif", function(x, value) {
    x@match_motif <- value
    validObject(x)
    x
})

#' @rdname CoRecMotif-class
setMethod("match_pvalue<-", "CoRecMotif", function(x, value) {
    x@match_pvalue <- value
    validObject(x)
    x
})

#' @rdname CoRecMotif-class
setMethod("match_cluster<-", "CoRecMotif", function(x, value) {
    x@match_cluster <- value
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
            is(x@match_motif, "universalmotif"),
            x@match_motif@altname,
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
        "   Probe set:      ", object@probe_set, "\n",
        "   PBM condition:  ", object@pbm_condition, "\n",
        "   Array ID:       ", object@array_id, "\n",
        "   Motif strength: ", round(object@motif_strength, 2), "\n",
        "   Rolling IC:     ", round(object@rolling_ic, 2), "\n\n"
    )
    print(round(object@ppm@motif, 2))
})

