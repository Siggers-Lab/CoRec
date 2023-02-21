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

summarize_corecmotifs <- function(corecmotifs) {
    # Convert each corecmotif object into a data frame
    corecmotif_df <-
        lapply(corecmotifs, as.data.frame) %>%

        # Combine all the data frames
        dplyr::bind_rows()

    # Return the dataframe of corecmotif information
    return(corecmotif_df)
}


set_seed_name <- function(corecmotif, seed_name) {
    assertthat::assert_that(assertthat::is.string(seed_name))

    corecmotif@seed_name <- seed_name

    return(corecmotif)
}

set_pbm_condition <- function(corecmotif, pbm_condition) {
    assertthat::assert_that(assertthat::is.string(pbm_condition))

    corecmotif@pbm_condition <- pbm_condition

    return(corecmotif)
}

set_array_id <- function(corecmotif, array_id) {
    assertthat::assert_that(assertthat::is.string(array_id))

    corecmotif@array_id <- array_id

    return(corecmotif)
}

set_seed_sequence <- function(corecmotif, seed_sequence) {
    assertthat::assert_that(assertthat::is.string(seed_sequence))

    corecmotif@seed_sequence <- seed_sequence

    return(corecmotif)
}

set_motif_name <- function(corecmotif, motif_name) {
    assertthat::assert_that(assertthat::is.string(motif_name))

    corecmotif@ppm@name <- motif_name

    return(corecmotif)
}

