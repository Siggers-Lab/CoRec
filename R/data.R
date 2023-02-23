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
