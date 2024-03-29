#' Plot CoRecMotifs and their matching reference motifs
#'
#' Creates a motif logo plot of a [CoRecMotif][CoRecMotif-class] and/or its
#' matching reference motif. There are several types of logos to choose from.
#'
#' @param corecmotif [CoRecMotif][CoRecMotif-class]. The
#'   [CoRecMotif][CoRecMotif-class] to plot.
#' @param corecmotif_logo_type `character(1)`. One of "delta_zscore", "ICM",
#'   "PWM", "PPM", or "none". The type of logo to plot for the
#'   [CoRecMotif][CoRecMotif-class]. (Default: "delta_zscore")
#' @param reference_logo_type `character(1)`. One of "ICM", "PWM", "PPM", or
#'   "none". The type of logo to plot for the matching reference motif.
#'   (Default: "ICM")
#' @param reverse_complement `logical(1)`. Should the CoRecMotif be reversed?
#'   (Default: False)
#' @param correct_orientation `logical(1)`. Should the reference motif be
#'   reversed if necessary to match the CoRecMotif's orientation? (Default:
#'   TRUE)
#' @param fade_corecmotif `logical(1)`. Should the
#'   [CoRecMotif][CoRecMotif-class] logo be faded? (Default: FALSE)
#' @param fade_reference `logical(1)`. Should the reference motif logo be faded?
#'   (Default: FALSE)
#' @param check_corecmotif `logical(1)`. Should `corecmotif` be checked for
#'   validity? Setting this to FALSE can increase speed, but if `corecmotif` is
#'   not a valid [CoRecMotif][CoRecMotif-class], it may produce uninformative
#'   error messages. (Default: TRUE)
#'
#' @return A `ggplot` object.
#'
#' @export
#'
#' @examples
#' print("FILL THIS IN")
plot_corecmotif <-
    function(
        corecmotif,
        corecmotif_logo_type = c("delta_zscore", "ICM", "PWM", "PPM", "none"),
        reference_logo_type = c("ICM", "PWM", "PPM", "none"),
        reverse_complement = FALSE,
        correct_orientation = TRUE,
        fade_corecmotif = FALSE,
        fade_reference = FALSE,
        check_corecmotif = TRUE
    ) {
    # Make sure all the arguments are the right type
    assertthat::assert_that(
        assertthat::is.flag(reverse_complement),
        assertthat::is.flag(correct_orientation),
        assertthat::is.flag(fade_corecmotif),
        assertthat::is.flag(fade_reference),
        assertthat::is.flag(check_corecmotif)
    )
    corecmotif_logo_type <- match.arg(corecmotif_logo_type)
    reference_logo_type <- match.arg(reference_logo_type)

    # Make sure corecmotif is a valid CoRecMotif
    if (check_corecmotif) {
        if (!methods::is(corecmotif, "CoRecMotif")) {
            stop(
                "corecmotif is not a CoRecMotif",
                call. = FALSE
            )
        }
        methods::validObject(corecmotif)
    }

    # Get the correct form of the CoRecMotif as a numeric matrix
    corecmotif_matrix <-
        switch(
            corecmotif_logo_type,
            "delta_zscore" = get_delta_zscore_motif(corecmotif),
            "ICM" = get_icm(corecmotif),
            "PWM" = get_pwm(corecmotif),
            "PPM" = get_ppm(corecmotif),
            "none" = NA
        )

    # Take the reverse complement if necessary
    if (reverse_complement && !is.na(corecmotif_matrix)) {
        corecmotif_matrix <-
            matrix(rev(corecmotif_matrix), nrow = 4, byrow = FALSE) %>%
            magrittr::set_rownames(c("A", "C", "G", "T")) %>%
            magrittr::set_colnames(1:ncol(.))
    }

    # Get the correct form of the reference motif as a numeric matrix
    reference_matrix <-
        switch(
            reference_logo_type,
            "ICM" = get_match_icm(corecmotif, correct_orientation),
            "PWM" = get_match_pwm(corecmotif, correct_orientation),
            "PPM" = get_match_ppm(corecmotif, correct_orientation),
            "none" = NA
        )

    # Take the reverse complement if necessary
    if (reverse_complement && correct_orientation && !is.na(reference_matrix)) {
        reference_matrix <-
            matrix(rev(reference_matrix), nrow = 4, byrow = FALSE) %>%
            magrittr::set_rownames(c("A", "C", "G", "T")) %>%
            magrittr::set_colnames(1:ncol(.))
    }

    # Set both plots to NA by default
    corecmotif_plot <- NA
    reference_plot <- NA

    # Plot the CoRecMotif if necessary
    if (is.matrix(corecmotif_matrix)) {
        corecmotif_plot <-
            plot_motif(
                corecmotif_matrix,
                logo_type = corecmotif_logo_type,
                outline_color = corecmotif_outline_color,
                fade_motif = fade_corecmotif
            ) +

            # Add a title
            ggplot2::ggtitle(
                paste0(
                    get_motif_name(corecmotif),
                    "\nMotif strength: ",
                    round(get_motif_strength(corecmotif), 3),
                    "\nRolling IC: ",
                    round(get_rolling_ic(corecmotif), 3)
                )
            )
    }

    # Plot the reference motif if necessary
    if (is.matrix(reference_matrix)) {
        reference_plot <-
            plot_motif(
                reference_matrix,
                logo_type = reference_logo_type,
                outline_color = reference_outline_color,
                fade_motif = fade_reference
            ) +

            # Add a title
            ggplot2::ggtitle(
                paste0(
                    get_match_altname(corecmotif),
                    "\nMatch p-value: ",
                    signif(get_match_pvalue(corecmotif), 3),
                    "\nMatch cluster: ",
                    get_match_cluster(corecmotif)
                )
            )
    }

    # Check if the motifs got plotted
    corecmotif_is_plot <- methods::is(corecmotif_plot, "gg")
    reference_is_plot <- methods::is(reference_plot, "gg")

    if (!corecmotif_is_plot && !reference_is_plot) {
        # If no plots were made, return NULL
        return()
    } else if (corecmotif_is_plot && !reference_is_plot) {
        # If there is no reference motif plot, return the CoRecMotif plot
        return(corecmotif_plot)
    } else if (!corecmotif_is_plot && reference_is_plot) {
        # If there is no CoRecMotif plot, return the reference motif plot
        return(reference_plot)
    } else {
        # If both plots were made, put them side by side and return
        combined_plot <-
            cowplot::plot_grid(corecmotif_plot, reference_plot, nrow = 1)

        return(combined_plot)
    }
}

# This is a private helper for plot_corecmotif()
plot_motif <- function(motif_matrix, logo_type, outline_color, fade_motif) {
    # Figure out what the y axis label should be based on the motif type
    y_label <-
        switch(
            logo_type,
            "delta_zscore" = expression(paste(Delta, "z-score", sep = "")),
            "PPM" = "Probability",
            "PWM" = "Weight",
            "ICM" = "Bits"
        )

    # Set the minimum and maximum values for the x axis
    xmin <- 0.5
    xmax <- ncol(motif_matrix) + 0.5

    # Make a small buffer so the tops of letters don't get cut off
    y_buffer = 0.03

    # Set the minimum and maximum values for the y axis
    if (logo_type == "ICM") {
        ymin <-0 - y_buffer
        ymax <- 2 + y_buffer
    } else if (logo_type == "PPM") {
        ymin <- 0 - y_buffer
        ymax <- 1 + y_buffer
    } else {
        # Sum all the negative numbers in each column and find the most extreme
        ymin <- min(sum_negatives(motif_matrix)) - y_buffer

        # Sum all the positive numbers in each column and find the most extreme
        ymax <- max(sum_positives(motif_matrix)) + y_buffer
    }

    # Start the motif logo plot
    motif_plot <-
        universalmotif::view_logo(
            motif_matrix,
            colour.scheme = logo_color_scheme,
            sort.positions = TRUE
        ) +

        # Add the y axis label
        ggplot2::ylab(y_label) +

        # Set the y axis limits
        ggplot2::ylim(c(ymin, ymax)) +

        # Add an outline to the positive section of the logo
        ggplot2::annotate(
            "rect",
            xmin = xmin,
            xmax = xmax,
            ymin = 0 - y_buffer,
            ymax = ymax,
            color = outline_color,
            fill = NA,
            linewidth = 1.5
        ) +

        # Set the formatting for the y axis text
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(),
            axis.title.y = ggplot2::element_text(angle = 90),
        )

    # Outline and fade the negative section of the logo (if it exists)
    if (logo_type %in% c("delta_zscore", "PWM")) {
        motif_plot <-
            motif_plot +

            # Add a semi-transparent box with an outline
            ggplot2::annotate(
                "rect",
                xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = 0,
                alpha = 0.65,
                colour = outline_color,
                fill = "white",
                linewidth = 1.5
            )
    }

    # Fade the whole logo if necessary
    if (fade_motif) {
        motif_plot <-
            motif_plot +

            # Add a semi-transparent box to the whole logo
            ggplot2::annotate(
                "rect",
                xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax,
                alpha = 0.65,
                colour = outline_color,
                fill = "white",
                linewidth = 1.5
            )
    }

    # Return the finished plot
    return(motif_plot)
}

