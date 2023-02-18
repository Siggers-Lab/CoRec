#' Make a list of corecmotif objects
#'
#' Creates a list of corecmotif objects using fluorescence data from a CoRec
#' experiment.
#'
#' This is a convenience function that calls
#' \code{\link{annotate_fluorescence_table}},
#' \code{\link{fluorescence_to_zscore_table}}, and
#' \code{\link{zscore_table_to_corecmotifs}} and optionally saves the outputs
#' using file names generated from user input. By default no output files are
#' created. To save output files, you must provide \code{output_directory},
#' \code{output_base_name}, or both. If \code{output_directory} is provided but
#' \code{output_base_name} is not, the files will be saved in the provided
#' directory using the base name "output". If \code{output_base_name} is
#' provided but \code{output_directory} is not, the files will be saved in the
#' current working directory using the provided base name. If \code{array_id} is
#' provided, it will be appended to \code{output_base_name}. The outputs of
#' \code{\link{annotate_fluorescence_table}},
#' \code{\link{fluorescence_to_zscore_table}}, and
#' \code{\link{zscore_table_to_corecmotifs}} will be saved with the suffixes
#' "_fluorescence.tsv", "_zscores.tsv", and "_all_corecmotifs.rds" respectively.
#'
#' @param fluorescence_file the path to the file containing the fluorescence
#'   data to load. See 'Details' of \code{\link{annotate_fluorescence_table}}
#'   for expected columns.
#' @param pbm_conditions a character vector specifying the PBM conditions (e.g.,
#'   cell type, treatment, and factor profiled) in the order they appear in
#'   \code{fluorescence_file}.
#' @param annotation_file the path to the file containing the probe annotations
#'   to use. See 'Details' of \code{\link{annotate_fluorescence_table}} for
#'   expected columns.
#' @param output_directory the directory where output files will be saved. No
#'   output files will be created unless output_directory or output_base_name is
#'   provided.
#' @param output_base_name the base name for output files. No output files will
#'   be created unless output_directory or output_base_name is provided.
#' @param array_id an optional (but recommended) tag specifying the particular
#'   array/experiment the fluorescence data is from.
#'
#' @return A list of \linkS4class{corecmotif} objects, one for each possible
#'   combination of the probe sets annotated in \code{annotation_file} and the
#'   PBM conditions listed in \code{pbm_conditions}.
#'
#' @export
#'
#' @seealso \link{annotate_fluorescence_table},
#'   \link{fluorescence_to_zscore_table}, \link{zscore_table_to_corecmotifs}
#' @examples
#' corecmotifs <-
#'     make_corecmotifs(
#'         fluorescence_file = "./example_data/hTF_v1_example_fluorescence.dat",
#'         pbm_conditions = c(
#'             "UT_SUDHL4_SMARCA4MIX",
#'             "UT_SUDHL4_HDAC1MIX",
#'              "UT_SUDHL4_SUZ12",
#'              "UT_SUDHL4_PRMT5"
#'         ),
#'         annotation_file = "example_data/hTF_v1_annotation.tsv",
#'         output_directory = "./example_data/example_output",
#'         output_base_name = "example",
#'         array_id = "v1_a6_run1"
#'     )
make_corecmotifs <-
    function(
        fluorescence_file,
        pbm_conditions,
        annotation_file,
        output_directory = NULL,
        output_base_name = NULL,
        array_id = NULL
    ) {
    # Do not save any output files by default
    fluorescence_output <- NULL
    zscore_output <- NULL
    corec_output <- NULL

    # If only the output directory or base name is provided, set the other
    if (is.null(output_directory) & !is.null(output_base_name)) {
        # If no output directory is provided, use current working directory
        output_directory <- getwd()
    } else if (!is.null(output_directory) & is.null(output_base_name)) {
        # If no output base name is provided, use "output"
        output_base_name <- "output"
    }

    # If both the output directory and base name exist, set output file names
    if (!is.null(output_directory) & !is.null(output_base_name)) {
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

        # Create names for the three output files
        fluorescence_output <- paste0(output_base_name, "_fluorescence.tsv")
        zscore_output <- paste0(output_base_name, "_zscores.tsv")
        corec_output <- paste0(output_base_name, "_all_corecmotifs.rds")
    }

    # Load and annotate the table of fluorescence values
    fluorescence_table <-
        annotate_fluorescence_table(
            fluorescence_file = fluorescence_file,
            pbm_conditions = pbm_conditions,
            annotation_file = annotation_file,
            output_file = fluorescence_output
        )

    # Convert the fluorescence values into condition-wise z-scores
    zscore_table <-
        fluorescence_to_zscore_table(
            fluorescence_table = fluorescence_table,
            fluorescence_columns = pbm_conditions,
            output_file = zscore_output
        )

    # Make corecmotif objects for all the seed/condition combos
    corecmotifs <-
        zscore_table_to_corecmotifs(
            zscore_table = zscore_table,
            zscore_columns = pbm_conditions,
            output_file = corec_output,
            array_id = array_id
        )

    # Return the list of corecmotifs
    return(corecmotifs)
}

# #' @param reference_motifs_file
# #' @param motif_strength_threshold
# #' @param rolling_ic_threshold
# #' @param comparison_method
# #' @param cluster_assignments_file
# #' @param pvalue_threshold
#
# motif_strength_threshold = 1,
# rolling_ic_threshold = 1.5,
# comparison_method = "ed",
# cluster_assignments_file = NULL,
# pvalue_threshold = 0.05
# reference_motifs_file,
#
# # Filter out corecmotifs with low scores
# filtered_corec_motifs <-
#     lapply(
#         corec_motifs,
#         function(corec_motif) {
#             # If the score is above the threshold, keep this motif
#             if (corec_motif@motif_strength > motif_strength_threshold &
#                 corec_motif@rolling_ic > rolling_ic_threshold) {
#                 return(corec_motif)
#             }
#         }
#     ) %>%
#
#     # Remove NULL elements (from corecmotifs below the score threshold)
#     plyr::compact()
#
# # Read in the table of cluster assignments (if provided)
# if (!is.null(cluster_assignments_file)) {
#     cluster_assignments <-
#         read.table(
#             cluster_assignments_file,
#             header = TRUE, sep = "\t"
#         )
# } else {
#     cluster_assignments <- NULL
# }
#
# # Compare the filtered corecmotifs to the library of reference TF motifs
# matched_corec_motifs <-
#     purrr::map(
#         filtered_corec_motifs,
#         identify_motif_match,
#         reference_motifs_file = reference_motifs_file,
#         cluster_assignments = cluster_assignments,
#         method = comparison_method
#     )
#
# # Save the list of matched corecmotifs as a single RDS file
# saveRDS(
#     matched_corec_motifs,
#     paste0(output_base_name, "_matched_corecmotifs.rds")
# )
#
# # Pull out the corecmotif PPMs that matched a reference motif well
# lapply(matched_corec_motifs, function(corec_motif) {
#     if (!(is.na(corec_motif@motif_match_pvalue)) &
#         corec_motif@motif_match_pvalue < pvalue_threshold) {
#         # Pull out the PPM
#         ppm <- corec_motif@ppm
#
#         # Add the name of the matched motif to the PPM name
#         ppm@name <-
#             paste0(ppm@name, "_", corec_motif@motif_match@altname)
#
#         # Return the updated PPM
#         return(ppm)
#     }
# }) %>%
#
#     # Remove NULL elements (from corecmotifs above the pvalue threshold)
#     plyr::compact() %>%
#
#     # Save a MEME format file of the matched corecmotifs
#     universalmotif::write_meme(
#         paste0(output_base_name, "_motifs.meme"),
#         overwrite = TRUE
#     )
#
# # Return
# return(matched_corec_motifs)


