# Seed probe z-score
#
# Finds the fluorescence value z-score of the seed probe of the z-score motif.
#
# @param zscore_motif A data frame representing a z-score motif, where the rows
#   are nucleotides and the columns are positions in the motif.
#
# @return The fluorescence value z-score of the seed probe of the z-score motif.
find_seed_zscore <- function(zscore_motif) {
    # The seed probe z-score shows up at every position of the z-score motif
    seed_probe_zscore <-
        Reduce(dplyr::intersect, as.list(zscore_motif))

    # Return the z-score
    return(seed_probe_zscore)
}


# Title
#
# @param ppm
# @param width
#
# @return
calculate_rolling_ic <- function(ppm, width = 5) {
    # Convert the PPM to an information content matrix
    icm <-
        universalmotif::convert_type(ppm, type = "ICM")

    #
    max_sliding_window_ic <-
        icm["motif"] %>%
        colSums() %>%
        zoo::rollmean(width) %>%
        max()

    names(max_sliding_window_ic) <- width

    return(max_sliding_window_ic)
}

