#' Title
#'
#' @param corec_motif
#' @param motif_score_type
#' @param motif_score_threshold
#'
#' @return
#' @export
#'
#' @examples
plot_corecmotif <-
    function(
        corec_motif,
        motif_logo_type = c("delta_zscore", "PPM", "PWM", "ICM"),
        motif_score_threshold = NA
    ) {
        # Make sure the selected motif logo type is a valid option
        motif_logo_type <-
            match.arg(motif_logo_type)

        # Get the correct form of the motif as a numeric matrix
        motif_matrix <-
            # Get the correct form of the motif
            switch(
                motif_logo_type,
                "delta_zscore" = corec_motif@delta_zscore_motif,
                "PPM" = corec_motif@ppm@motif,
                "PWM" =
                    universalmotif::convert_type(corec_motif@ppm, "PWM")@motif,
                "ICM" =
                    universalmotif::convert_type(corec_motif@ppm, "ICM")@motif
            ) %>%

            # Convert to a numeric matrix
            as.matrix()

        # Figure out the width of the motif
        motif_width <-
            ncol(corec_motif@zscore_motif)

        # Start the motif logo plot
        motif_plot <-
            ggseqlogo::ggseqlogo(
                motif_matrix,
                method = "custom",
                seq_type = "dna"
            )

        # Get the range of the y axis
        motif_plot_yrange <-
            ggplot2::ggplot_build(motif_plot)$layout$panel_params[[1]]$y.range

        # Format the width and height of the axes without printing any messages
        suppressMessages(
            motif_plot <-
                motif_plot +

                # Increase the plot width by decreasing the padding on the edges
                # ggseqlogo already makes an x scale, so this prints a message
                ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +

                # Increase the plot height and include 1 decimal place in labels
                ggplot2::scale_y_continuous(
                    expand = c(0.01, 0.01),
                    labels = function(x) {sprintf("%.1f", x)}
                )
        )

        # Add semi-transparent boxes to fade out parts of the logo
        motif_plot <-
            motif_plot +

            # Add a semi-transparent box over the negative section of the logo
            ggplot2::annotate(
                'rect',
                xmin = 0.5,
                xmax = motif_width + 0.5,
                ymin = motif_plot_yrange[1],
                ymax = 0.0,
                alpha = 0.65,
                col='white',
                fill='white'
            ) +

            # Add another box over the whole logo if the motif score is too low
            if (!is.na(motif_score_threshold) &
                slot(corec_motif, "motif_score") < motif_score_threshold
            ) {
                ggplot2::annotate(
                    'rect',
                    xmin = 0.5,
                    xmax = motif_width + 0.5,
                    ymin = motif_plot_yrange[1],
                    ymax= motif_plot_yrange[2],
                    alpha = 0.65,
                    col='white',
                    fill='white'
                )
            }

        # Add a blue outline to the negative section of the logo if necessary
        if (motif_logo_type %in% c("delta_zscore", "PWM")) {
            motif_plot <-
                motif_plot +

                ggplot2::annotate(
                    'rect',
                    xmin = 0.5,
                    xmax = motif_width + 0.5,
                    ymin = motif_plot_yrange[1],
                    ymax = 0,
                    color = "#abbed1",
                    fill = NA,
                    size = 1.5
                )
        }

        # Figure out what the y axis label should be based on the motif type
        y_label <-
            switch(
                motif_logo_type,
                "delta_zscore" = expression(paste(Delta, "z-score", sep="")),
                "PPM" = "probability",
                "PWM" = "weight",
                "ICM" = "bits"
            )

        # Add remaining plot elements
        motif_plot <-
            motif_plot +

            # Add a blue outline to the positive section of the logo
            ggplot2::annotate(
                'rect',
                xmin = 0.5,
                xmax = motif_width + 0.5,
                ymin = 0,
                ymax = motif_plot_yrange[2],
                color = "#abbed1",
                fill = NA,
                size = 1.5
            ) +

            # Add a title with the PBM condition, seed name, and motif score
            ggplot2::ggtitle(
                paste(
                    paste("PBM condition:", corec_motif@pbm_condition),
                    paste("Seed:", corec_motif@seed_name),
                    paste0(
                        slot(corec_motif, "motif_score_type"),
                        ": ",
                        round(slot(corec_motif, "motif_score"), digits = 3)
                    ),
                    motif_logo_type,
                    sep="\n"
                )
            ) +

            # Add the y axis label
            ggplot2::ylab(y_label) +

            # Set the formatting for the axis text and the plot title
            ggplot2::theme(
                axis.text.x = ggplot2::element_blank(),
                axis.title.x = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_text(size = 10, face = "plain"),
                axis.title.y = ggplot2::element_text(size = 10, face = "plain"),
                plot.title = ggplot2::element_text(size = 10)
            )

        # Return the finished plot
        return(motif_plot)
    }

