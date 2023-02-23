# Paste together the output directory, output base name, and array ID
update_output_base_name <-
    function(
        output_directory = NULL,
        output_base_name = NULL,
        array_id = NULL
    ) {
    if (is.null(output_directory) && is.null(output_base_name)) {
        # Return NULL if neither the output directory or base name is provided
        return(NULL)
    } else if (is.null(output_directory) && !is.null(output_base_name)) {
        # If no output directory is provided, use current working directory
        output_directory <- getwd()
    } else if (!is.null(output_directory) && is.null(output_base_name)) {
        # If no output base name is provided, use "output"
        output_base_name <- "output"
    }

    # Create the output directory if it doesn't already exist
    if (!dir.exists(output_directory)) {
        dir.create(output_directory)
    }

    # Add the output directory path to the base file name for output files
    output_base_name <- paste(output_directory, output_base_name, sep = "/")

    # Add the array ID to the base file name if it's provided
    if (!is.null(array_id)) {
        output_base_name <- paste(output_base_name, array_id, sep = "_")
    }

    # Return the updated output base name
    return(output_base_name)
}

# Try to save a TSV or RDS file but catch errors and give a warning instead
try_catch_save_output <- function(x, output_file, file_type = c("tsv", "rds")) {
    # Make sure the file type is either tsv or rds
    file_type <- match.arg(file_type)

    # Save the file only if an output file path is provided
    if (!is.null(output_file)) {
        tryCatch(
            {
                if (file_type == "tsv") {
                    # Try to save a TSV file
                    suppressWarnings(
                        write.table(
                            x,
                            output_file,
                            quote = FALSE,
                            sep = "\t",
                            row.names = FALSE,
                            col.names = TRUE
                        )
                    )
                } else if (file_type == "rds") {
                    # Try to save an RDS file
                    suppressWarnings(
                        saveRDS(
                            x,
                            output_file
                        )
                    )
                }
            },
            # If it fails, skip the output saving step with a warning
            error = function(e) {
                warning(
                    "Could not write to output file '",
                    output_file,
                    "'\nSkipping output file creation...",
                    call. = FALSE
                )
            }
        )
    }
}

# Figure out the seed probe z-score from a z-score motif
# This is separate from get_seed_zscore because this is used during creation of
#   a new CoRecMotif, meaning you don't have the whole CoRecMotif object yet
find_seed_zscore <- function(zscore_motif) {
    seed_zscore <-
        # The seed probe z-score shows up at every position of the z-score motif
        which(table(zscore_motif) >= ncol(zscore_motif)) %>%

        # The names of the frequency table are the z-scores
        names() %>%

        # Convert to numeric
        as.numeric()

    # Return the z-score
    return(seed_zscore)
}

# Calculate the beta parameter to use when converting z-score motifs to PPMs
calculate_beta <- function(zscore_motif) {
    # Find the seed probe z-score
    seed_zscore <- find_seed_zscore(zscore_motif)

    # Calculate beta
    beta <- 4 - (0.5 * seed_zscore)

    # Restrict beta to a range of 1 to 4
    beta <- max(min(4, beta), 1)

    # Return beta
    return(beta)
}

# Calculate the motif strength (median of top 15% of probes)
calculate_strength <- function(zscore_motif) {
    # Get a sorted list of all the probe z-scores
    z_scores <-
        sort(unique(unlist(zscore_motif)), decreasing = TRUE)

    # Figure out how many probes to average (rounded up to the nearest integer)
    num_probes <- ceiling(length(z_scores) * (15 / 100))

    # Find the median z-score of the highest num_probes probes
    median_zscore <- median(z_scores[1:num_probes])

    # Return the average z-score
    return(median_zscore)
}

# Convert a z-score motif to a PPM
zscore_to_ppm <- function(zscore_motif, beta, name = "motif") {
    # Transform the z-scores using the beta parameter
    ppm <-
        # Multiply each z-score by beta and then take the exponential
        exp(beta * zscore_motif) %>%

        # Normalize by dividing each value in each column by the column-wise sum
        apply(2, function(col) {col / sum(col)}) %>%

        # Convert to a universalmotif object
        universalmotif::create_motif(name = name)

    # Return the PPM
    return(ppm)
}

# Calculate the rolling information content over a window of length 5
calculate_rolling_ic <- function(ppm) {
    # Convert the PPM to an information content matrix
    icm <-
        universalmotif::convert_type(ppm, type = "ICM")

    #
    max_sliding_window_ic <-
        icm@motif %>%
        colSums() %>%
        zoo::rollmean(5) %>%
        max()

    return(max_sliding_window_ic)
}

