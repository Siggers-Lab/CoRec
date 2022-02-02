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
        seed_zscore_threshold = 1,
        rolling_ic_threshold = 1.5,
        outline_color = "#abbed1"
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

        # Fade the whole logo if either motif score is below its threshold
        fade_motif <-
            (!is.na(seed_zscore_threshold) &
             corec_motif@seed_zscore < seed_zscore_threshold) |
            (!is.na(rolling_ic_threshold) &
             corec_motif@rolling_ic < rolling_ic_threshold)

        # Make the motif logo plot
        motif_plot <-
            plot_motif_helper(
                motif_matrix,
                motif_logo_type,
                outline_color,
                fade_motif
            )

        # Add a title with the PBM condition, seed name, and motif score
        motif_plot <-
            motif_plot +

            ggplot2::ggtitle(
                paste(
                    paste("PBM condition:", corec_motif@pbm_condition),
                    paste("Seed:", corec_motif@seed_name),
                    paste0(
                        "Seed z-score: ",
                        round(corec_motif@seed_zscore, 3)
                    ),
                    paste0(
                        "Rolling IC: ",
                        round(corec_motif@rolling_ic, 3)
                    ),
                    motif_logo_type,
                    sep="\n"
                )
            )

        # Return the finished plot
        return(motif_plot)
    }


plot_reference_motif <-
    function(
        corec_motif,
        motif_logo_type = c("ICM", "PPM", "PWM"),
        outline_color = "#added1"
    ) {
        # Make sure the selected motif logo type is a valid option
        motif_logo_type <-
            match.arg(motif_logo_type)

        # Get the correct form of the reference motif as a numeric matrix
        motif_matrix <-
            # Get the correct form of the motif
            switch(
                motif_logo_type,
                "PPM" = corec_motif@motif_match@motif,
                "PWM" =
                    universalmotif::convert_type(
                        corec_motif@motif_match,
                        "PWM"
                    )@motif,
                "ICM" =
                    universalmotif::convert_type(
                        corec_motif@motif_match,
                        "ICM"
                    )@motif
            ) %>%

            # Convert to a numeric matrix
            as.matrix()

        # Make the motif logo plot
        motif_plot <-
            plot_motif_helper(
                motif_matrix,
                motif_logo_type,
                outline_color
            )

        # Add a title with the name of the reference motif and the logo type
        motif_plot <-
            motif_plot +

            ggplot2::ggtitle(
                paste(
                    "Best match reference motif",
                    paste(
                        corec_motif@motif_match@name,
                        corec_motif@motif_match@altname,
                        sep = "_"
                    ),
                    paste(
                        "Match p-value:",
                        signif(corec_motif@motif_match_pvalue, 3)
                    ),
                    motif_logo_type,
                    "",
                    sep="\n"
                )
            )

        # Return the finished plot
        return(motif_plot)
    }


plot_corec_reference_motifs <-
    function(
        corec_motif,
        motif_logo_type = c("delta_zscore", "PPM", "PWM", "ICM"),
        seed_zscore_threshold = 1,
        rolling_ic_threshold = 1.5,
        reference_logo_type = c("ICM", "PPM", "PWM")
    ) {
        # Make sure the selected motif logo types are valid options
        motif_logo_type <-
            match.arg(motif_logo_type)
        reference_logo_type <-
            match.arg(reference_logo_type)

        # Plot the corecmotif
        corec_motif_plot <-
            plot_corecmotif(
                corec_motif,
                motif_logo_type,
                seed_zscore_threshold,
                rolling_ic_threshold,
                outline_color = "#abbed1"
            )

        # Plot the reference motif
        reference_motif_plot <-
            plot_reference_motif(
                corec_motif,
                reference_logo_type,
                outline_color = "#added1"
            )

        # Combine the corecmotif and the reference motif into one plot
        combined_plot <-
            cowplot::plot_grid(
                corec_motif_plot,
                reference_motif_plot,
                nrow = 1
            )

        # Return the final plot
        return(combined_plot)
    }


plot_motif_helper <-
    function(
        motif,
        motif_logo_type = c("delta_zscore", "PPM", "PWM", "ICM"),
        outline_color = "#abbed1",
        fade_motif = FALSE
    ) {
        # Make sure the selected motif logo type is a valid option
        motif_logo_type <-
            match.arg(motif_logo_type)

        # Figure out the width of the motif
        motif_width <-
            ncol(motif)

        # Start the motif logo plot
        # Don't print the dumb "`guides(<scale> = FALSE)` is deprecated" warning
        # It's not my fault. It's in the ggseqlogo source code. :'(
        suppressWarnings(
            motif_plot <-
                ggseqlogo::ggseqlogo(
                    motif,
                    method = "custom",
                    seq_type = "dna"
                )
        )

        # Make the y-range go from 0 to 2 if it's an ICM plot
        if (motif_logo_type == "ICM") {
            motif_plot <-
                motif_plot +

                ggplot2::ylim(c(0,2))
        }

        # Get the range of the y axis
        motif_plot_yrange <-
            ggplot2::ggplot_build(motif_plot)$layout$panel_params[[1]]$y.range

        # Format the width and height of the axes
        # Don't print the stupid "Scale for 'x' is already present" message
        # Once again, not my fault. It's the ggseqlogo source code. Again. :'(
        suppressMessages(
            motif_plot <-
                motif_plot +

                # Increase the plot width by decreasing the padding on the edges
                # ggseqlogo already makes an x scale, so this prints a message
                ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +

                # Increase the plot height and include 1 decimal place in labels
                # This prints a message if it's an ICM plot (that is my fault)
                ggplot2::scale_y_continuous(
                    expand = c(0.01, 0.01),
                    labels = function(x) {sprintf("%.1f", x)}
                )
        )

        # Fade the whole logo if necessary
        if (fade_motif) {
            motif_plot <-
                motif_plot +

                ggplot2::annotate(
                    'rect',
                    xmin = 0.5,
                    xmax = motif_width + 0.5,
                    ymin = motif_plot_yrange[1],
                    ymax = motif_plot_yrange[2],
                    alpha = 0.65,
                    col='white',
                    fill='white'
                )
        }

        # Outline and fade the negative section of the logo (if it exists)
        if (motif_logo_type %in% c("delta_zscore", "PWM")) {
            motif_plot <-
                motif_plot +

                # Add a semi-transparent box
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

                # Add an outline
                ggplot2::annotate(
                    'rect',
                    xmin = 0.5,
                    xmax = motif_width + 0.5,
                    ymin = motif_plot_yrange[1],
                    ymax = 0,
                    color = outline_color,
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

            # Add an outline to the positive section of the logo
            ggplot2::annotate(
                'rect',
                xmin = 0.5,
                xmax = motif_width + 0.5,
                ymin = 0,
                ymax = motif_plot_yrange[2],
                color = outline_color,
                fill = NA,
                size = 1.5
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

