get_seed_name <- function(corecmotif) {
    return(corecmotif@seed_name)
}

get_pbm_condition <- function(corecmotif) {
    return(corecmotif@pbm_condition)
}

get_array_id <- function(corecmotif) {
    return(corecmotif@array_id)
}

get_beta <- function(corecmotif) {
    return(corecmotif@beta)
}

get_zscore_motif <- function(corecmotif) {
    return(corecmotif@zscore_motif)
}

get_rolling_ic <- function(corecmotif) {
    return(corecmotif@rolling_ic)
}

get_motif_strength <- function(corecmotif) {
    return(corecmotif@motif_strength)
}

get_seed_sequence <- function(corecmotif) {
    return(corecmotif@seed_sequence)
}

# Seed probe z-score
#
# Finds the fluorescence value z-score of the seed probe of the z-score motif.
#
# @param zscore_motif A data frame representing a z-score motif, where the rows
#   are nucleotides and the columns are positions in the motif.
#
# @return The fluorescence value z-score of the seed probe of the z-score motif.
get_seed_zscore <- function(corecmotif) {
    # Get the z-score motif
    zscore_motif <- get_zscore_motif(corecmotif)

    seed_zscore <- find_seed_zscore(zscore_motif)

    # Return the z-score
    return(seed_zscore)
}

get_ppm <- function(corecmotif) {
    return(corecmotif@ppm)
}

get_motif_name <- function(corecmotif) {
    return(corecmotif@ppm@name)
}

get_match_motif <- function(corecmotif) {
    return(corecmotif@match)
}

get_match_name <- function(corecmotif) {
    return(corecmotif@match@name)
}

get_match_altname <- function(corecmotif) {
    return(corecmotif@match@altname)
}

get_match_pvalue <- function(corecmotif) {
    return(corecmotif@match_pvalue)
}

get_match_cluster <- function(corecmotif) {
    return(corecmotif@match_cluster)
}

#' Make a data frame summarizing a list of CoRecMotifs
#'
#' Creates a data frame representation of a list of
#' [CoRecMotifs][CoRecMotif-class].
#'
#' @param corecmotifs the list of [CoRecMotifs][CoRecMotif-class] to summarize.
#'
#' @return A data frame with information about a list of
#'   [CoRecMotifs][CoRecMotif-class].
#' @export
#'
#' @examples
#' print("FILL THIS IN")
summarize_corecmotifs <- function(corecmotifs) {
    # Convert each corecmotif object into a data frame
    corecmotif_df <-
        lapply(corecmotifs, as.data.frame) %>%

        # Combine all the data frames
        dplyr::bind_rows()

    # Return the dataframe of corecmotif information
    return(corecmotif_df)
}

