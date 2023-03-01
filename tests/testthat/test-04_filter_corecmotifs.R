test_that("bad arguments to character parameters are handled correctly", {
    # All fail when given a value that isn't a character vector
    expect_error(
        filter_corecmotifs(
            example_corecmotifs,
            probe_set = list()
        ),
        "probe_set is not a character vector"
    )

    expect_error(
        filter_corecmotifs(
            example_corecmotifs,
            pbm_condition = FALSE
        ),
        "pbm_condition is not a character vector"
    )

    expect_error(
        filter_corecmotifs(
            example_corecmotifs,
            seed_sequence = NA
        ),
        "seed_sequence is not a character vector"
    )

    expect_error(
        filter_corecmotifs(
            example_corecmotifs,
            motif_name = matrix()
        ),
        "motif_name is not a character vector"
    )

    expect_error(
        filter_corecmotifs(
            example_corecmotifs,
            match_name = data.frame()
        ),
        "match_name is not a character vector"
    )

    expect_error(
        filter_corecmotifs(
            example_corecmotifs,
            match_altname = c(5, 3, 1)
        ),
        "match_altname is not a character vector"
    )

    expect_error(
        filter_corecmotifs(
            example_corecmotifs,
            match_cluster = 10
        ),
        "match_cluster is not a character vector"
    )
})

test_that("bad arguments to numeric parameters are handled correctly", {
    # All fail when given a value that isn't numeric
    expect_error(
        filter_corecmotifs(
            example_corecmotifs,
            motif_strength = list()
        ),
        "motif_strength is not a number"
    )

    expect_error(
        filter_corecmotifs(
            example_corecmotifs,
            rolling_ic = FALSE
        ),
        "rolling_ic is not a number"
    )

    expect_error(
        filter_corecmotifs(
            example_corecmotifs,
            match_pvalue = NA
        ),
        "match_pvalue is not a number"
    )

    # All fail when given a numeric vector with length > 1
    expect_error(
        filter_corecmotifs(
            example_corecmotifs,
            motif_strength = c(1, 2, 3)
        ),
        "motif_strength is not a number"
    )

    expect_error(
        filter_corecmotifs(
            example_corecmotifs,
            rolling_ic = c(3, 2, 1)
        ),
        "rolling_ic is not a number"
    )

    expect_error(
        filter_corecmotifs(
            example_corecmotifs,
            match_pvalue = c(NA_real_, 0)
        ),
        "match_pvalue is not a number"
    )
})

test_that("filtering by character parameters works", {
    # Filter by probe_set ------------------------------------------------------
    probe_set_filter <-
        filter_corecmotifs(
            example_corecmotifs,
            probe_set = c("MA0079.3_SP1", "MA1418.1_IRF3")
        )

    expect_equal(probe_set_filter, example_corecmotifs[c(1:4, 17:24, 37:40)])

    # Filter by pbm_condition --------------------------------------------------
    pbm_condition_filter <-
        filter_corecmotifs(
            example_corecmotifs,
            pbm_condition = "UT_SUDHL4_HDAC_mix"
        )

    expect_equal(
        pbm_condition_filter,
        example_corecmotifs[c(2, 6, 10, 14, 18, 22, 26, 30, 34, 38)]
    )

    # Filter by seed_sequence --------------------------------------------------
    seed_sequence_filter <-
        filter_corecmotifs(
            example_corecmotifs,
            seed_sequence =
                "GCTCCAAGATGGCGGCGTCGATCGAGCACGCAGATCGTCTTGATTCGCTTGACGCTGCTG"
        )

    expect_equal(seed_sequence_filter, example_corecmotifs[c(5:8, 25:28)])

    # Filter by motif_name -----------------------------------------------------
    motif_name_filter <-
        filter_corecmotifs(
            example_corecmotifs,
            motif_name = c(
                "MA0095.2_YY1_UT_SUDHL4_JMJD2A_v1_a11_run1",
                "MA1418.1_IRF3_UT_SUDHL4_SWISNF_mix_v1_a11_run1"
            )
        )

    expect_equal(motif_name_filter, example_corecmotifs[c(8, 17)])

    # Filter by match_name -----------------------------------------------------
    match_name_filter <-
        filter_corecmotifs(
            example_matched_corecmotifs,
            match_name = "MA0755.1"
        )

    expect_equal(match_name_filter, example_matched_corecmotifs[5:8])

    # Filter by match_altname --------------------------------------------------
    match_altname_filter <-
        filter_corecmotifs(
            example_matched_corecmotifs,
            match_altname = c("MA0517.1.STAT1::STAT2", "MA1651.1.ZFP42")
        )

    expect_equal(match_altname_filter, example_matched_corecmotifs[c(1, 9)])

    # Filter by match_cluster --------------------------------------------------
    match_cluster_filter <-
        filter_corecmotifs(
            example_matched_corecmotifs,
            match_cluster = "MEF"
        )

    expect_equal(match_cluster_filter, example_matched_corecmotifs[3:4])
})

test_that("filtering by numeric parameters works", {
    # Filter by motif_strength -------------------------------------------------
    motif_strength_filter <-
        filter_corecmotifs(
            example_corecmotifs,
            motif_strength = 2
        )

    expect_equal(
        motif_strength_filter,
        example_corecmotifs[c(1, 4, 8, 11, 14, 19, 24, 28, 31, 33:34, 39:40)]
    )

    # Filter by rolling_ic -----------------------------------------------------
    rolling_ic_filter <-
        filter_corecmotifs(
            example_corecmotifs,
            rolling_ic = 1.6
        )

    expect_equal(
        rolling_ic_filter,
        example_corecmotifs[c(1, 8, 13:14, 18, 24, 27:28, 31, 33:34, 39:40)]
    )

    # Filter by match_pvalue ---------------------------------------------------
    match_pvalue_filter <-
        filter_corecmotifs(
            example_matched_corecmotifs,
            match_pvalue = 1e-5
        )

    expect_equal(
        match_pvalue_filter,
        example_matched_corecmotifs[c(1, 3:4, 9)]
    )
})

test_that("filtering by multiple parameters works", {
    # Filter by probe_set and pbm_condition ------------------------------------
    combined_filter_1 <-
        filter_corecmotifs(
            example_corecmotifs,
            probe_set = c("MA0079.3_SP1", "MA1418.1_IRF3"),
            pbm_condition = c("UT_SUDHL4_HDAC_mix", "UT_SUDHL4_PRMT5")
        )

    expect_equal(
        combined_filter_1,
        example_corecmotifs[c(2:3, 18:19, 22:23, 38:39)]
    )

    # Filter by pbm_condition and motif_strength -------------------------------
    combined_filter_2 <-
        filter_corecmotifs(
            example_corecmotifs,
            pbm_condition = c("UT_SUDHL4_HDAC_mix", "UT_SUDHL4_PRMT5"),
            motif_strength = 1
        )

    expect_equal(
        combined_filter_2,
        example_corecmotifs[c(3, 11, 14, 18:19, 23, 31, 34, 38:39)]
    )

    # Filter by motif_strength and rolling_ic ----------------------------------
    combined_filter_3 <-
        filter_corecmotifs(
            example_corecmotifs,
            motif_strength = 1,
            rolling_ic = 1.5
        )

    expect_equal(
        combined_filter_3,
        example_corecmotifs[c(1, 8, 11, 13:14, 18, 24, 28, 31, 33:34, 39:40)]
    )
})

