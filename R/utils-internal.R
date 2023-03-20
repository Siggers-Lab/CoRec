# Set useful constants ---------------------------------------------------------

# Set the colors to use for nucleotides when plotting motif logos
logo_color_scheme <- c(
    "A" = "#109648",
    "C" = "#255C99",
    "G" = "#F7B32B",
    "T" = "#D62839"
)

# Set the colors to use for motif logo outlines
corecmotif_outline_color <- "#abbed1"
reference_outline_color <- "#added1"

# Define a bunch of dummy variables to get rid of R CMD check notes
motif <- motif_1 <- motif_2 <- NULL
probe_id <- probe_sequence <- probe_set <- NULL
snv_position <- snv_nucleotide <- NULL
cluster <- n_conditions <- distance <- NULL
pbm_condition <- pbm_conditions <- zscore <- NULL
array_id <-  motif_strength <- rolling_ic <- NULL
match_cluster <- match_motif <- match_pvalue <- NULL

# Define useful functions ------------------------------------------------------

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

# Make sure a data frame has the expected column names
check_colnames <- function(x, expected_columns) {
    # Make sure all the expected columns are present
    if (!all(expected_columns %in% colnames(x))) {
        stop(
            deparse(substitute(x)),
            " is missing one or more expected columns\n",
            "Expected columns: ",
            paste(expected_columns, collapse = ", "),
            call. = FALSE
        )
    }

    # Remove any extra columns
    x <- dplyr::select(x, dplyr::all_of(expected_columns))

    # Also get rid of duplicate rows if there are any
    x <- dplyr::distinct(x)

    # Return the data frame
    return(x)
}

# Make sure an object is a valid z-score motif
check_valid_zscore_motif <- function(zscore_motif) {
    # Make sure it's a numeric matrix
    assertthat::assert_that(
        is.matrix(zscore_motif) && is.numeric(zscore_motif)
    )

    # Make sure there are four rows
    if (nrow(zscore_motif) != 4) {
        stop(
            "zscore_motif must have four rows\n",
            "Rows must be in the order 'A', 'C', 'G', 'T'",
            call. = FALSE
        )
    }

    # Check if the row names are the nucleotides and fix them if not
    if (!identical(rownames(zscore_motif), c("A", "C", "G", "T"))) {
        warning(
            "zscore_motif row names are not 'A', 'C', 'G', 'T'\n",
            "Make sure rows are in the order 'A', 'C', 'G', 'T'\n",
            "Changing row names.",
            call. = FALSE
        )
        rownames(zscore_motif) <- c("A", "C", "G", "T")
    }

    # Make sure the column names are "1", "2", "3", etc.
    colnames(zscore_motif) <- as.character(1:ncol(zscore_motif))

    return(zscore_motif)
}

check_corecmotif_list <- function(corecmotifs) {
    # If it's not a list or an individual CoRecMotif, give an error
    if (!is.list(corecmotifs) && !methods::is(corecmotifs, "CoRecMotif")) {
        stop(
            "corecmotifs is not a list of CoRecMotifs or coercible to one",
            call. = FALSE
        )
    }

    # If it's an individual valid CoRecMotif, return it in a list
    if (!is.list(corecmotifs) && methods::is(corecmotifs, "CoRecMotif")) {
        methods::validObject(corecmotifs)
        return(list(corecmotifs))
    }

    # Check if each element of corecmotifs is a CoRecMotif
    are_corecmotifs <-
        vapply(corecmotifs, methods::is, class2 = "CoRecMotif", logical(1))

    # If any elements of corecmotifs are not CoRecMotifs, give an error
    if (!all(are_corecmotifs)) {
        stop(
            "not all elements of corecmotifs are CoRecMotifs",
            call. = FALSE
        )
    } else if (length(corecmotifs) < 1) {
        stop(
            "corecmotifs is an empty list",
            call. = FALSE
        )
    }

    # Check if each element of corecmotifs is a VALID CoRecMotif
    vapply(corecmotifs, methods::validObject, logical(1))
    return(corecmotifs)
}

# Try to save a TSV or RDS file but catch errors and give a warning instead
try_catch_save_output <- function(x, output_file, file_type = c("tsv", "rds")) {
    # Make sure the file type is either tsv or rds
    file_type <- match.arg(file_type)

    # Make sure the output_file argument is a string
    assertthat::assert_that(
        assertthat::is.string(output_file) || is.null(output_file)
    )

    # Save the file only if an output file path is provided
    if (!is.null(output_file)) {
        tryCatch(
            {
                if (file_type == "tsv") {
                    # Try to save a TSV file
                    suppressWarnings(
                        utils::write.table(
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

# Create an array ID if one isn't provided
create_array_id <- function() {
    # Paste 8 random digits after "random_id_"
    array_id <-
        paste(
            "random_id",
            paste(sample(0:9, 8, replace = TRUE), collapse = ""),
            sep = "_"
        )

    # Return the ID
    return(array_id)
}

# Figure out the seed probe z-score from a z-score motif
# This is separate from the seed_zscore method because this is used during
# creation of a new CoRecMotif, meaning you don't have the whole object yet
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
calculate_beta <- function(motif_strength) {
    # Calculate beta
    beta <- 4 - (0.5 * motif_strength)

    # Restrict beta to a range of 1 to 4
    beta <- max(min(4, beta), 1)

    # Return beta
    return(beta)
}

# Calculate the motif strength (93rd percentile of z-scores)
calculate_strength <- function(zscore_motif) {
    motif_strength <-
        zscore_motif %>%
        sort() %>%
        unique() %>%
        stats::quantile(0.93, names = FALSE)

    # Return the motif strength
    return(motif_strength)
}

# Convert a z-score motif to a PPM
zscore_to_universalmotif <- function(zscore_motif, beta, name = "motif") {
    # Transform the z-scores using the beta parameter
    motif <-
        # Multiply each z-score by beta and then take the exponential
        exp(beta * zscore_motif) %>%

        # Normalize by dividing each value in each column by the column-wise sum
        apply(2, function(col) {col / sum(col)}) %>%

        # Convert to a universalmotif object
        universalmotif::create_motif(name = name)

    # Return the PPM
    return(motif)
}

# Calculate the max mean rolling information content over a window of length 5
calculate_rolling_ic <- function(motif) {
    # Convert the PPM to an information content (IC) matrix
    icm <- universalmotif::convert_type(motif, type = "ICM")

    # Figure out the total IC at each position
    ic_per_position <- universalmotif::colSums(icm)

    # Get the mean IC over each window of length 5
    max_rolling_ic <-
        vapply(
            1:(length(ic_per_position) - 5),
            function(i) {
                mean(ic_per_position[i:(i + 4)])
            },
            numeric(1)
        ) %>%

        # Find the window with the highest mean IC
        max()

    # Return the max rolling IC
    return(max_rolling_ic)
}

# Sum all the negative numbers in each column of a matrix
sum_negatives <- function(motif_matrix) {
    apply(motif_matrix, 2, function(column) {
        sum(column[column < 0])
    })
}

# Sum all the positive numbers in each column of a matrix
sum_positives <- function(motif_matrix) {
    apply(motif_matrix, 2, function(column) {
        sum(column[column > 0])
    })
}

