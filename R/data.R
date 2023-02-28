#' hTF Array annotation (v1)
#'
#' A table containing probe IDs and metadata for Version 1 of the hTF Array.
#'
#' @format A data frame with 17,434 rows and 6 columns:
#' \describe{
#'   \item{probe_id}{The probe ID.}
#'   \item{probe_type}{The probe type, either "MOTIF" or "BACKGROUND".}
#'   \item{probe_sequence}{The probe sequence.}
#'   \item{probe_set}{The name of the probe set this probe is part of.}
#'   \item{snv_position}{The position of the SNV in the sequence.}
#'   \item{snv_nucleotide}{The nucleotide at the SNV position.}
#' }
"hTF_v1_annotation"

#' Motif cluster assignments
#'
#' A table assigning reference motifs to clusters.
#'
#' @format A data frame with 946 rows and 2 columns:
#' \describe{
#'   \item{motif}{The name of a reference motif.}
#'   \item{cluster}{The name of the cluster this reference motif is part of.}
#' }
"motif_clusters"

#' Example fluorescence and z-score data
#'
#' Tables of example CoRec fluorescence data and z-score data for 5 probe sets
#' in 4 conditions.
#'
#' @format A data frame with 533 rows and 5 or 10 columns:
#'
#'   ## example_fluorescence_table
#'
#'   The un-annotated fluorescence table contains 5 columns: `probe_id`, and 4
#'   columns of raw fluorescence data named `UT_SUDHL4_SWISNF_mix`,
#'   `UT_SUDHL4_HDAC_mix`, `UT_SUDHL4_PRMT5`, and `UT_SUDHL4_JMJD2A`.
#'
#'   ## example_annotated_fluorescence_table
#'
#'   The annotated fluorescence table contains the 6 probe information columns
#'   described in [hTF_v1_annotation] and the 4 raw fluorescence data columns
#'   described above.
#'
#'   ## example_zscore_table
#'
#'   The z-score table contains the same columns as the annotated fluorescence
#'   table after converting the raw fluorescence data into z-scores.
#'
#' @rdname example_fluorescence_tables
"example_fluorescence_table"

#' @rdname example_fluorescence_tables
#' @format NULL
"example_annotated_fluorescence_table"

#' @rdname example_fluorescence_tables
#' @format NULL
"example_zscore_table"

#' Example CoRecMotifs
#'
#' Lists of example [CoRecMotifs][CoRecMotif-class].
#'
#' @format A list of 10 or 40 [CoRecMotifs][CoRecMotif-class]
#'
#'   ## example_corecmotifs
#'
#'   A list of 40 [CoRecMotifs][CoRecMotif-class] that have not been matched to
#'   reference motifs.
#'
#'   ## example_matched_corecmotifs
#'
#'   A list of 10 [CoRecMotifs][CoRecMotif-class] that have been matched to
#'   reference motifs.
#'
#' @rdname example_corecmotifs
"example_corecmotifs"

#' @rdname example_corecmotifs
#' @format NULL
"example_matched_corecmotifs"
