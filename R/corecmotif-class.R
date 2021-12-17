# Import the universalmotif class to use it in a class union
#' @importClassesFrom universalmotif universalmotif

# Define a class union that can contain either a universalmotif object or NULL
setClassUnion("universalmotif_or_NULL", c("universalmotif", "NULL"))

#' Title
#'
#' @slot seed_name character.
#' @slot pbm_condition character.
#' @slot beta numeric.
#' @slot zscore_motif data.frame.
#' @slot delta_zscore_motif data.frame.
#' @slot ppm universalmotif.
#' @slot motif_score_type character.
#' @slot motif_score numeric.
#' @slot seed_probe_sequence character.
#' @slot motif_match universalmotif_or_NULL.
#' @slot motif_match_score_type
#' @slot motif_match_score
#' @slot motif_match_pvalue
#'
#' @return
#' @export
#'
#' @examples
setClass(
    # Name the class CoRecMotif
    "corecmotif",

    # Define the names and types of the slots the class should have
    slots = list(
        seed_name = "character",
        pbm_condition = "character",
        beta = "numeric",
        zscore_motif = "data.frame",
        delta_zscore_motif = "data.frame",
        ppm = "universalmotif",
        motif_score_type = "character",
        motif_score = "numeric",
        motif_strength = "numeric",
        seed_probe_sequence = "character",
        motif_match = "universalmotif_or_NULL",
        motif_match_score_type = "character",
        motif_match_score = "numeric",
        motif_match_pvalue = "numeric",
        motif_match_qvalue = "numeric",
        motif_cluster_match = "character"
    ),

    # Provide a default example object
    prototype = list(
        seed_name = NA_character_,
        pbm_condition = NA_character_,
        beta = NA_real_,
        zscore_motif = data.frame(),
        delta_zscore_motif = data.frame(),
        ppm = universalmotif::create_motif("ACGT"),
        motif_score_type = NA_character_,
        motif_score = NA_real_,
        motif_strength = NA_real_,
        seed_probe_sequence = NA_character_,
        motif_match = NULL,
        motif_match_score_type = NA_character_,
        motif_match_score = NA_real_,
        motif_match_pvalue = NA_real_,
        motif_match_qvalue = NA_real_,
        motif_cluster_match = NA_character_
    )
) %>%

    # Do not print the message in universalmotif's constructor
    suppressMessages()



#' Title
#'
#' @param seed_name
#' @param pbm_condition
#' @param beta
#' @param zscore_motif
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
corecmotif <-
    function(
        seed_name,
        pbm_condition,
        zscore_motif,
        beta_method = "linear",
        ...
    ) {
        delta_zscore_motif <-
            zscore_to_delta_zscore(zscore_motif)

        motif_name <-
            paste(seed_name, pbm_condition, sep = "_")

        beta <-
            calculate_beta(zscore_motif, method = beta_method)

        ppm <-
            zscore_to_ppm(zscore_motif, beta, motif_name)

        new(
            "corecmotif",
            seed_name = seed_name,
            pbm_condition = pbm_condition,
            beta = beta,
            zscore_motif = zscore_motif,
            delta_zscore_motif = delta_zscore_motif,
            ppm = ppm,
            ...
        )
    }

