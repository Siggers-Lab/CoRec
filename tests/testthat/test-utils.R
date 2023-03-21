test_that("summarize_corecmotifs() works", {
    expected_df <-
        data.frame(
            "probe_set" = c(
                example_corecmotifs[[2]]@probe_set,
                example_corecmotifs[[1]]@probe_set,
                example_matched_corecmotifs[[3]]@probe_set
            ),
            "pbm_condition" = c(
                example_corecmotifs[[2]]@pbm_condition,
                example_corecmotifs[[1]]@pbm_condition,
                example_matched_corecmotifs[[3]]@pbm_condition
            ),
            "array_id" = c(
                example_corecmotifs[[2]]@array_id,
                example_corecmotifs[[1]]@array_id,
                example_matched_corecmotifs[[3]]@array_id
            ),
            "list_index" = c(
                2,
                1,
                3
            ),
            "motif_strength" = c(
                example_corecmotifs[[2]]@motif_strength,
                example_corecmotifs[[1]]@motif_strength,
                example_matched_corecmotifs[[3]]@motif_strength
            ),
            "rolling_ic" = c(
                example_corecmotifs[[2]]@rolling_ic,
                example_corecmotifs[[1]]@rolling_ic,
                example_matched_corecmotifs[[3]]@rolling_ic
            ),
            "seed_sequence" = c(
                example_corecmotifs[[2]]@seed_sequence,
                example_corecmotifs[[1]]@seed_sequence,
                example_matched_corecmotifs[[3]]@seed_sequence
            ),
            "match_motif" = c(
                NA,
                NA,
                example_matched_corecmotifs[[3]]@match_motif["altname"]
            ),
            "match_pvalue" = c(
                NA_real_,
                NA_real_,
                example_matched_corecmotifs[[3]]@match_pvalue
            ),
            "match_cluster" = c(
                NA_character_,
                NA_character_,
                example_matched_corecmotifs[[3]]@match_cluster
            ),
            "best_match_cluster" = c(
                NA_character_,
                NA_character_,
                example_matched_corecmotifs[[3]]@match_cluster
            )
        ) %>%

        # Convert to a tibble
        tibble::as_tibble()

    expect_equal(
        summarize_corecmotifs(
            list(
                example_corecmotifs[[1]],
                example_corecmotifs[[2]],
                example_matched_corecmotifs[[3]]
            )
        ),
        expected_df
    )
})

test_that("bad arguments for cluster_assignments are handled correctly", {
    # Fails when given a value that isn't a data frame
    expect_error(
        update_cluster_match(
            example_matched_corecmotifs[[1]],
            cluster_assignments =
                c("motif_1" = "cluster_1", "motif_2" = "cluster_2")
        ),
        "cluster_assignments is not a data frame"
    )

    # Fails when given a data frame with the wrong column names
    expect_error(
        update_cluster_match(
            example_matched_corecmotifs[[1]],
            cluster_assignments =
                data.frame(
                    "V1" = c("motif_1", "motif_2"),
                    "V2" = c("cluster_1", "cluster_2")
                )
        ),
        "cluster_assignments is missing one or more expected columns"
    )

    # Gives a warning when cluster_assignments is missing the motif name
    expect_warning(
        update_cluster_match(
            example_matched_corecmotifs[[1]],
            cluster_assignments =
                data.frame(
                    "motif" = c("motif_1", "motif_2"),
                    "cluster" = c("cluster_1", "cluster_2")
                )
        ),
        "Motif match altname not in cluster assignments table; "
    )
})

test_that("update_cluster_match() works", {
    motif_1 <- example_matched_corecmotifs[[1]]
    motif_1@match_motif <- get_match_motif(example_matched_corecmotifs[[5]])

    expect_equal(
        get_match_cluster(
            update_cluster_match(motif_1, cluster_assignments = motif_clusters)
        ),
        get_match_cluster(example_matched_corecmotifs[[5]])
    )

    expect_equal(
        get_match_cluster(
            update_cluster_match(motif_1)
        ),
        NA_character_
    )

    expect_equal(
        get_match_cluster(
            update_cluster_match(
                example_corecmotifs[[1]], cluster_assignments = motif_clusters
            )
        ),
        NA_character_
    )
})
