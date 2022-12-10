#' Title
#'
#' @param output_directory
#' @param fluorescence_file
#' @param pbm_conditions
#' @param annotation_file
#' @param reference_motifs_file
#' @param array_id
#' @param motif_strength_threshold
#' @param rolling_ic_threshold
#' @param comparison_method
#' @param cluster_assignments_file
#' @param pvalue_threshold
#'
#' @return
#' @export
#'
#' @examples
run_full_analysis <-
    function(
        output_directory,
        fluorescence_file,
        pbm_conditions,
        annotation_file,
        reference_motifs_file,
        array_id = NA,
        motif_strength_threshold = 1,
        rolling_ic_threshold = 1.5,
        comparison_method = "ed",
        cluster_assignments_file = NULL,
        pvalue_threshold = 0.05
    ) {
        # Load and annotate the table of fluorescence values
        fluorescence_table <-
            make_fluorescence_table(
                fluorescence_file,
                pbm_conditions,
                annotation_file,
                array_id
            )

        # Update pbm_conditions based on the column names of fluorescence_table
        pbm_conditions <-
            fluorescence_table %>%

            # Select the columns whose names contain PBM condition names
            dplyr::select(dplyr::contains(pbm_conditions)) %>%

            # Keep just the column names
            names()

        # Create the output directory if it doesn't already exist
        if (!dir.exists(output_directory)) {
            dir.create(output_directory)
        }

        # Create a base file name for output files
        output_base_name <-
            paste0(
                output_directory,
                "/output"
            )

        # Add the array ID to the base file name if it's provided
        if (!is.na(array_id)) {
            output_base_name <-
                paste0(
                    output_base_name,
                    "_",
                    array_id
                )
        }

        # Save the annotated fluorescence matrix
        write.table(
            fluorescence_table,
            paste0(output_base_name, "_fluorescence.dat"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE
        )

        # Convert the fluorescence values into condition-wise z-scores
        zscore_matrix <-
            fluorescence_to_zscore_matrix(
                fluorescence_table,
                pbm_conditions
            )

        # Save the z-score matrix
        write.table(
            zscore_matrix,
            paste0(output_base_name, "_zscores.dat"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE
        )

        # Extract a vector of all the non-background probe seed names
        seed_names <-
            zscore_matrix %>%

            # Keep only one row (the seed probe row) for each TF probe set
            dplyr::filter(SNV_pos_offset == 0) %>%

            # Get just the column containing the seed names
            dplyr::pull(seed_names)

        # Get a list of corecmotif objects for all the seed/condition combos
        corec_motifs <-
            zscore_matrix_to_motifs(
                zscore_matrix,
                seed_names,
                pbm_conditions
            )

        # Create a directory to save the corecmotifs if it doesn't already exist
        if (!dir.exists(paste0(output_directory, "/corecmotifs_all"))) {
            dir.create(paste0(output_directory, "/corecmotifs_all"))
        }

        # Save all the corecmotifs as individual RDS files
        lapply(corec_motifs, function(corec_motif) {
            file_name <-
                paste0(
                    output_directory,
                    "/corecmotifs_all/",
                    corec_motif@ppm@name,
                    ".rds"
                )

            saveRDS(corec_motif, file_name)
        })

        # Save the list of all corecmotifs as a single RDS file
        saveRDS(
            corec_motifs,
            paste0(output_base_name, "_all_corecmotifs.rds")
        )

        # Filter out corecmotifs with low scores
        filtered_corec_motifs <-
            lapply(
                corec_motifs,
                function(corec_motif) {
                    # If the score is above the threshold, keep this motif
                    if (corec_motif@motif_strength > motif_strength_threshold &
                        corec_motif@rolling_ic > rolling_ic_threshold) {
                        return(corec_motif)
                    }
                }
            ) %>%

            # Remove NULL elements (from corecmotifs below the score threshold)
            plyr::compact()

        # Read in the table of cluster assignments (if provided)
        if (!is.null(cluster_assignments_file)) {
            cluster_assignments <-
                read.table(
                    cluster_assignments_file,
                    header = TRUE, sep = "\t"
                )
        } else {
            cluster_assignments <- NULL
        }

        # Compare the filtered corecmotifs to the library of reference TF motifs
        matched_corec_motifs <-
            purrr::map(
                filtered_corec_motifs,
                identify_motif_match,
                reference_motifs_file = reference_motifs_file,
                cluster_assignments = cluster_assignments,
                method = comparison_method
            )

        # Create a directory to save the corecmotifs with matched TF motifs
        if (!dir.exists(paste0(output_directory, "/corecmotifs_matched"))) {
            dir.create(paste0(output_directory, "/corecmotifs_matched"))
        }

        # Save the matched corecmotifs as individual RDS files
        lapply(matched_corec_motifs, function(corec_motif) {
            file_name <-
                paste0(
                    output_directory,
                    "/corecmotifs_matched/",
                    corec_motif@ppm@name,
                    ".rds"
                )

            saveRDS(corec_motif, file_name)
        })

        # Save the list of matched corecmotifs as a single RDS file
        saveRDS(
            matched_corec_motifs,
            paste0(output_base_name, "_matched_corecmotifs.rds")
        )

        # Pull out the corecmotif PPMs that matched a reference motif well
        lapply(matched_corec_motifs, function(corec_motif) {
            if (!(is.na(corec_motif@motif_match_pvalue)) &
                corec_motif@motif_match_pvalue < pvalue_threshold) {
                # Pull out the PPM
                ppm <- corec_motif@ppm

                # Add the name of the matched motif to the PPM name
                ppm@name <-
                    paste0(ppm@name, "_", corec_motif@motif_match@altname)

                # Return the updated PPM
                return(ppm)
            }
        }) %>%

            # Remove NULL elements (from corecmotifs above the pvalue threshold)
            plyr::compact() %>%

            # Save a MEME format file of the matched corecmotifs
            universalmotif::write_meme(
                paste0(output_base_name, "_motifs.meme")
            )

        # Return
        return(matched_corec_motifs)
    }

