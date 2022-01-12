extract_dataframe <- function(corec_motifs) {
    # Extract the relevant information from each corecmotif object
    corecmotif_list <- lapply(corec_motifs, function(corec_motif) {
        return(
            list(
                seed_name = corec_motif@seed_name,
                pbm_condition = corec_motif@pbm_condition,
                seed_zscore = corec_motif@seed_zscore,
                motif_strength = corec_motif@motif_strength,
                motif_match = ifelse(
                    is.null(corec_motif@motif_match),
                    NA,
                    corec_motif@motif_match@altname
                ),
                motif_match_pvalue = corec_motif@motif_match_pvalue,
                motif_cluster_match = corec_motif@motif_cluster_match
            )
        )
    })

    # Convert the list of corecmotif information into a dataframe
    corecmotif_dataframe <- dplyr::bind_rows(corecmotif_list)

    # Add a column with the seed TF name without the MA* identifier
    corecmotif_dataframe <-
        corecmotif_dataframe %>%

        # Remove everything up to and including the first underscore
        dplyr::mutate(seed_TF = gsub(".*_", "", seed_name)) %>%

        # Move the new column right after the seed_name column
        dplyr::relocate(seed_TF, .after = seed_name)

    # Return the dataframe of corecmotif information
    return(corecmotif_dataframe)
}

combine_replicates <- function(input_file, output_directory, output_base_name) {
    # Load the input file of replicate information
    replicates_dataframe <-
        read.table(input_file, header = FALSE, sep = "\t")

    # Read in the RDS files of lists of corecmotif objects
    corec_motifs <- lapply(unique(replicates_dataframe$V1), function(rds) {
        readRDS(rds)
    })

    # Rename each list of corecmotif objects to its corresponding RDS file name
    names(corec_motifs) <- unique(replicates_dataframe$V1)

    # Extract the dataframes of relevant info for each list of corecmotifs
    corec_motifs_dataframes <-
        lapply(corec_motifs, function(corec_motif_list) {
            extract_dataframe(corec_motif_list)
        })

    # Pull out the relevant rows for each replicate from the relevant dataframe
    corec_motifs_dataframes_subset <-
        lapply(1:nrow(replicates_dataframe), function(index) {
            # Get the relevant dataframe
            relevant_rows <-
                corec_motifs_dataframes[[replicates_dataframe$V1[index]]] %>%

                # Pull out the relevant rows
                dplyr::filter(pbm_condition == replicates_dataframe$V2[index])

            # Replace the values in the pbm_condition column with the group name
            relevant_rows$pbm_condition <- replicates_dataframe$V3[index]

            # Return the dataframe of relevant rows
            return(relevant_rows)
        })

    # Combine all the separate dataframes into one, averaging replicates
    combined_dataframe <-

        # Put all the rows from each dataframe into one dataframe
        dplyr::bind_rows(corec_motifs_dataframes_subset) %>%

        # Change the name of the PBM condition column to group
        dplyr::rename(group = pbm_condition) %>%

        # Group by seed name (and seed TF) and group
        dplyr::group_by(seed_name, seed_TF, group) %>%

        # Remove any groups that don't have 2 replicates to average
        dplyr::filter(dplyr::n() == 2) %>%

        # Average the seed z-score and motif strength for each seed across reps
        dplyr::summarise(
            seed_zscore = mean(seed_zscore),
            motif_strength = mean(motif_strength),
            motif_match_1 = motif_match[1],
            motif_match_pvalue_1 = motif_match_pvalue[1],
            motif_cluster_match_1 = motif_cluster_match[1],
            motif_match_2 = motif_match[2],
            motif_match_pvalue_2 = motif_match_pvalue[2],
            motif_cluster_match_2 = motif_cluster_match[2]
        )

    # Save the combined, averaged dataframe in long format
    write.table(
        combined_dataframe,
        paste0(output_directory, "/", output_base_name, "_long.tsv"),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE
    )

    # Convert the combined dataframe to wide format with the seed z-score
    combined_dataframe_seed <-
        combined_dataframe %>%

        # Remove the motif strength and match information columns
        dplyr::select(-motif_strength, -dplyr::contains("match")) %>%

        # Convert to wide format with the seed z-score as the values
        tidyr::pivot_wider(names_from = group, values_from = seed_zscore)

    # Save the wide format seed z-score matrix
    write.table(
        combined_dataframe_seed,
        paste0(
            output_directory, "/", output_base_name, "_wide_seed_zscore.tsv"
        ),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE
    )

    # Convert the combined dataframe to wide format with the motif strength
    combined_dataframe_strength <-
        combined_dataframe %>%

        # Remove the seed z-score and match information columns
        dplyr::select(-seed_zscore, -dplyr::contains("match")) %>%

        # Convert to wide format with the motif strength as the values
        tidyr::pivot_wider(names_from = group, values_from = motif_strength)

    # Save the wide format motif strength matrix
    write.table(
        combined_dataframe_strength,
        paste0(
            output_directory, "/", output_base_name, "_wide_motif_strength.tsv"
        ),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE
    )

    # Return the long format combined, averaged dataframe
    return(combined_dataframe)
}

