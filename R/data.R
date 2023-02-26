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

#' Example fluorescence data
#'
#' Tables of example CoRec fluorescence data for 5 probe sets in 4 conditions.
#'
#' @format A data frame with 533 rows and 5 columns:
#' \describe{
#'   \item{probe_id}{The probe ID.}
#'   \item{UT_SUDHL4_SWISNF_mix}{The fluorescence values for the
#'     "UT_SUDHL4_SWISNF_mix" PBM condition.}
#'   \item{UT_SUDHL4_HDAC_mix}{The fluorescence values for the
#'     "UT_SUDHL4_HDAC_mix" PBM condition.}
#'   \item{UT_SUDHL4_PRMT5}{The fluorescence values for the
#'     "UT_SUDHL4_PRMT5" PBM condition.}
#'   \item{UT_SUDHL4_JMJD2A}{The fluorescence values for the
#'     "UT_SUDHL4_SWISNF_mix" PBM condition.}
#' }
#'
#' @rdname example_fluorescence_tables
"example_fluorescence_table_1"

#' @rdname example_fluorescence_tables
#' @format NULL
"example_fluorescence_table_2"

