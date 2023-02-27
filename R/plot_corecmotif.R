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
        reference_logo_type = c("ICM", "PWM", "PPM", "none")
    ) {
    # Make sure the selected motif logo types are valid options
    corecmotif_logo_type <- match.arg(corecmotif_logo_type)
    reference_logo_type <- match.arg(reference_logo_type)

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

    # Get the correct form of the reference motif as a numeric matrix
    reference_matrix <-
        switch(
            reference_logo_type,
            "ICM" = get_match_icm(corecmotif),
            "PWM" = get_match_pwm(corecmotif),
            "PPM" = get_match_ppm(corecmotif),
            "none" = NA
        )

    # Set both plots to NA by default
    corecmotif_plot <- NA
    reference_plot <- NA

    # Plot the CoRecMotif if necessary
    if (is.matrix(corecmotif_matrix)) {
        corecmotif_plot <-
            plot_motif(
                corecmotif_matrix,
                logo_type = corecmotif_logo_type,
                outline_color = corecmotif_outline_color
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
                outline_color = reference_outline_color
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
plot_motif <- function(motif_matrix, logo_type, outline_color) {
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

    # Set the minimum and maximum values for the y axis
    if (logo_type == "ICM") {
        ymin <- -0.01
        ymax <- 2.01
    } else if (logo_type == "PPM") {
        ymin <- -0.01
        ymax <- 1.01
    } else {
        # Sum all the negative numbers in each column and find the most extreme
        ymin <- min(sum_negatives(motif_matrix)) - 0.01

        # Sum all the positive numbers in each column and find the most extreme
        ymax <- max(sum_positives(motif_matrix)) + 0.01
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
            'rect',
            xmin = xmin,
            xmax = xmax,
            ymin = 0,
            ymax = ymax,
            color = outline_color,
            fill = NA,
            size = 1.5
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
                'rect',
                xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = 0,
                alpha = 0.65,
                colour = outline_color,
                fill = "white",
                size = 1.5
            )
    }

    # Return the finished plot
    return(motif_plot)
}

