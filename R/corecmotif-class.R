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
#' @slot seed_zscore numeric.
#' @slot rolling_ic numeric.
#' @slot motif_strength numeric.
#' @slot seed_probe_sequence character.
#' @slot motif_match universalmotif_or_NULL.
#' @slot motif_match_method character.
#' @slot motif_match_pvalue numeric.
#' @slot motif_match_qvalue numeric.
#' @slot motif_cluster_match character.
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
        seed_zscore = "numeric",
        rolling_ic = "numeric",
        motif_strength = "numeric",
        seed_probe_sequence = "character",
        motif_match = "universalmotif_or_NULL",
        motif_match_method = "character",
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
        seed_zscore = NA_real_,
        rolling_ic = NA_real_,
        motif_strength = NA_real_,
        seed_probe_sequence = NA_character_,
        motif_match = NULL,
        motif_match_method = NA_character_,
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
        ic_window_width = 5,
        top_n_percent = 15,
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

        seed_zscore <-
            find_seed_zscore(zscore_motif)

        rolling_ic <-
            calculate_rolling_ic(ppm, ic_window_width)

        motif_strength <-
            calculate_strength(zscore_motif, top_n_percent)

        new(
            "corecmotif",
            seed_name = seed_name,
            pbm_condition = pbm_condition,
            beta = beta,
            zscore_motif = zscore_motif,
            delta_zscore_motif = delta_zscore_motif,
            ppm = ppm,
            seed_zscore = seed_zscore,
            rolling_ic = rolling_ic,
            motif_strength = motif_strength,
            ...
        )
    }

