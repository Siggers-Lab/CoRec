#' Title
#'
#' Identifies the reference motif that is the best match to the corecmotif and
#' calculates a p-value for that match.
#'
#' @param corec_motif
#' @param reference_motifs
#' @param method
#'
#' @return
#' @export
#'
#' @examples
identify_motif_match <-
    function(
        corec_motif,
        reference_motifs_file,
        cluster_assignments = NULL,
        min_overlap = 5,
        method = c(
            "ed",
            "allr",
            "pearson",
            "sandelin",
            "kullback",
            "blic1",
            "blic5",
            "llr1",
            "llr5"
        ),
        meme_path = "/share/pkg.7/meme/5.3.3/install/bin/"
    ) {
        # Make sure the selected comparison method is a valid option
        method <-
            match.arg(method)

        # Compare the corecmotif to the full library of reference motifs
        motif_comparison <-
            memes::runTomTom(
                corec_motif@ppm,
                database = reference_motifs_file,
                thresh = 1,
                evalue = FALSE,
                min_overlap = min_overlap,
                dist = method,
                meme_path = meme_path
            )

        # If no matches were found, return the original corecmotif object
        if (is.na(motif_comparison$best_match_motif)) {
            return(corec_motif)
        }

        # Update the motif matching slots of the original corecmotif object
        corec_motif@motif_match <-
            motif_comparison$best_match_motif[[1]]

        corec_motif@motif_match_method <-
            method

        corec_motif@motif_match_pvalue <-
            as.numeric(motif_comparison$best_match_pval)

        corec_motif@motif_match_qvalue <-
            as.numeric(motif_comparison$best_match_qval)

        # If there are no cluster assignments, return the updated corecmotif
        if (is.null(cluster_assignments)) {
            return(corec_motif)
        }

        # Figure out the cluster with the lowest average p-value-based distance
        best_cluster <-
            motif_comparison$tomtom[[1]] %>%

            # Convert the p-value into a distance metric of sorts
            dplyr::mutate(pval_distance = pvalue_to_distance(match_pval)) %>%

            # Merge with the dataframe of cluster information
            dplyr::full_join(
                cluster_assignments,
                by = c("match_altname" = "motif")
            ) %>%

            # Group by cluster
            dplyr::group_by(cluster) %>%

            # Find the average p-value-based distance for each cluster
            dplyr::summarise(mean_distance = mean(pval_distance)) %>%

            # Find the cluster with the smallest average distance
            dplyr::slice_min(order_by = mean_distance) %>%

            # Pull out just the cluster name
            dplyr::pull(cluster)

        # Update the cluster match slot of the original corecmotif object
        corec_motif@motif_cluster_match <-
            as.character(best_cluster)

        # Return the updated corecmotif
        return(corec_motif)
    }


cluster_motifs <-
    function(
        motifs_file,
        min_overlap = 5,
        method = c(
            "ed",
            "allr",
            "pearson",
            "sandelin",
            "kullback",
            "blic1",
            "blic5",
            "llr1",
            "llr5"
        ),
        meme_path = "/share/pkg.7/meme/5.3.3/install/bin/"
    ) {
        # Make sure the selected comparison method is a valid option
        method <-
            match.arg(method)

        # Compare all of the motifs to each other
        motif_comparison <-
            memes::runTomTom(
                motifs_file,
                database = motifs_file,
                thresh = 1,
                evalue = FALSE,
                min_overlap = min_overlap,
                dist = method,
                meme_path = meme_path
            )

        # Pull out the motif name, match name, and match p-value information
        motif_distance_matrix <-
            lapply(1:nrow(motif_comparison), function(row_index) {
                # Get the motif name
                motif_name <- motif_comparison$altname[row_index]

                # Get the motif match name and p-value distance metric columns
                motif_dataframe <-
                    motif_comparison$tomtom[[row_index]] %>%

                    # Convert the p-value into a distance metric of sorts
                    dplyr::mutate(
                        pval_distance = pvalue_to_distance(match_pval)
                    ) %>%

                    # Keep only the match name and the p-value distance metric
                    dplyr::select(
                        match_altname,
                        pval_distance
                    )

                # Add the motif name to the dataframe of match information
                motif_dataframe$motif_name <- motif_name

                # Return the dataframe
                return(motif_dataframe)
            }) %>%

            # Combine all the dataframes into one
            dplyr::bind_rows() %>%

            # Convert the dataframe into wide format
            tidyr::pivot_wider(
                names_from = match_altname,
                values_from = pval_distance
            ) %>%

            # Convert the motif names into rownames
            tibble::column_to_rownames("motif_name")

        # Convert the distance matrix into a dist object
        motif_distance_matrix <-
            dist(motif_distance_matrix)

        # Perform hierarchical clustering on the motifs
        clustered_motifs <-
            hclust(motif_distance_matrix, method = "complete")

        # Return the hierarchically clustered motifs
        return(clustered_motifs)
    }

pvalue_to_distance <- function(pvalue) {
    # Take the -log10 of the p-value
    log_pval = -log10(pvalue)

    # Cap the -log10 of the p-value at 15
    log_pval_capped = pmin(log_pval, 15)

    # Convert to a distance instead of a similarity
    pval_distance = 15 - log_pval_capped

    # Return the p-value-based distance metric
    return(pval_distance)
}

plot_motif_cluster <- function(motif_names, motif_list, pdf_filename = NA) {
    # Pull out the names of the motifs in motif_list
    motif_list_names <-
        lapply(motif_list, function(motif) {
            return(motif@altname)
        }) %>%
        unlist()

    # Pull out the universalmotif objects corresponding to motif_names
    motifs <- motif_list[which(motif_list_names %in% motif_names)]

    # Make a plot of each individual universalmotif
    motif_plots <-
        lapply(motifs, function(motif) {
            ggseqlogo::ggseqlogo(motif@motif) +
                ggplot2::ggtitle(motif@altname)
        })

    # Figure out how many columns and rows to have
    number_cols <- ceiling(sqrt(length(motif_names)))
    number_rows <- ceiling(length(motif_names) / number_cols)

    # Combine the individual motif plots in a grid
    combined_plot <-
        cowplot::plot_grid(
            plotlist = motif_plots,
            ncol = number_cols,
            nrow = number_rows
        )

    # Save the combined plot as a pdf if a file name is given
    if (!is.na(pdf_filename)) {
        ggplot2::ggsave(
            pdf_filename,
            plot = combined_plot,
            width = number_cols * 3,
            height = number_rows * 1.5,
            units = "in"
        )
    }

    # Return the combined plot
    return(combined_plot)
}

plot_all_motif_clusters <-
    function(cluster_assignments, motif_list, output_directory) {
        # Make a PDF file of the motifs in each cluster
        invisible(lapply(unique(cluster_assignments), function(cluster) {
            # Find the names of the motifs in this cluster
            motif_names <- names(which(cluster_assignments == cluster))

            # Generate a file name for this cluster
            pdf_filename <- paste0(
                output_directory,
                "/cluster_",
                cluster,
                "_motifs.pdf"
            )

            # Make and save the plot
            plot_motif_cluster(motif_names, motif_list, pdf_filename)
        }))
    }


plot_dendrogram <-
    function(
        clustered_motifs,
        cluster_assignments
    ) {
        # Convert the cluster assignments vector into a two column data frame
        cluster_assignments_df <-
            data.frame(
                label = names(cluster_assignments),
                cluster = factor(cluster_assignments)
            )

        # Convert to a format ggplot can use
        clustered_motifs <-
            ggdendro::dendro_data(clustered_motifs, type = "rectangle")

        # Add the cluster numbers to the labels dataframe
        clustered_motifs$labels <-
            dplyr::left_join(
                clustered_motifs$labels,
                cluster_assignments_df,
                by = "label"
            )

        # Add the clusters to the leaf node segments and renumber them
        leaf_node_segments <-
            clustered_motifs$segments %>%

            # Keep only leaf nodes
            dplyr::filter(yend == 0) %>%

            # Merge with the labels dataframe by the x coordinate
            dplyr::left_join(
                clustered_motifs$labels, by = "x"
            ) %>%

            # Renumber the clusters in order from left to right
            dplyr::rowwise() %>%
            dplyr::mutate(
                renumbered_cluster = which(unique(.$cluster) == cluster)
            )

        # Convert the renumbered clusters column to a factor
        leaf_node_segments$renumbered_cluster <-
            as.factor(leaf_node_segments$renumbered_cluster)

        # Add the renumbered cluster to the labels dataframe
        clustered_motifs$labels <-
            dplyr::left_join(
                clustered_motifs$labels,
                dplyr::select(leaf_node_segments, cluster, renumbered_cluster),
                by = "cluster"
            )

        # Start the dendrogram plot
        ggplot2::ggplot() +

            # Add the line segments
            ggplot2::geom_segment(
                data = ggdendro::segment(clustered_motifs),
                ggplot2::aes(x = x, y = y, xend = xend, yend = yend)
            ) +

            # Color the leaf node line segments
            ggplot2::geom_segment(
                data = leaf_node_segments,
                ggplot2::aes(
                    x = x,
                    y = y.x,
                    xend = xend,
                    yend = yend,
                    color = renumbered_cluster
                )
            ) +

            # Add the labels
            ggplot2::geom_text(
                data = ggdendro::label(clustered_motifs),
                ggplot2::aes(
                    x,
                    y,
                    label = label,
                    hjust= 1,
                    vjust = 0.5,
                    color = renumbered_cluster,
                    angle = 90
                ),
                size=3
            ) +

            # Set the colorscheme
            ggplot2::scale_color_manual(
                values = rep(RColorBrewer::brewer.pal(5, "Set1"), times = 100)
                ) +

            # Set the theme
            ggplot2::theme_bw() +

            # Remove unnecessary elements
            ggplot2::theme(
                panel.grid = ggplot2::element_blank(),
                rect = ggplot2::element_blank(),
                legend.position = "none",
                axis.title = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank()
            )
    }
