test_that("bad arguments for zscore_table are handled correctly", {
    # Check zscore_table_to_corecmotifs() --------------------------------------

    # Fails when given a value that isn't a data frame
    expect_error(
        zscore_table_to_corecmotifs(
            4,
            zscore_columns = c(
                "UT_SUDHL4_SWISNF_mix", "UT_SUDHL4_HDAC_mix"
            )
        ),
        "zscore_table is not a data frame"
    )

    # Fails when given a data frame with the wrong column names
    expect_error(
        zscore_table_to_corecmotifs(
            data.frame("probe_id" = c("1", "2", "3"), "column_2" = c(4, 5, 6)),
            zscore_columns = c(
                "UT_SUDHL4_SWISNF_mix", "UT_SUDHL4_HDAC_mix"
            )
        ),
        "zscore_table is missing one or more expected columns"
    )

    # Check make_zscore_motif() ------------------------------------------------

    # Fails when given a value that isn't a data frame
    expect_error(
        make_zscore_motif(
            4,
            probe_set_name = "MA0079.3_SP1",
            pbm_condition = "UT_SUDHL4_SWISNF_mix"
        ),
        "zscore_table is not a data frame"
    )

    # Fails when given a data frame with the wrong column names
    expect_error(
        make_zscore_motif(
            data.frame("probe_id" = c("1", "2", "3"), "column_2" = c(4, 5, 6)),
            probe_set_name = "MA0079.3_SP1",
            pbm_condition = "UT_SUDHL4_SWISNF_mix"
        ),
        "zscore_table is missing one or more expected columns"
    )

    # Fails when given a data frame with the wrong number of seed probe rows
    expect_error(
        make_zscore_motif(
            dplyr::filter(example_zscore_table, snv_position > 0),
            probe_set_name = "MA0079.3_SP1",
            pbm_condition = "UT_SUDHL4_SWISNF_mix"
        ),
        "Expected 1 seed probe for the probe set"
    )
})

test_that("bad arguments for zscore_columns are handled correctly", {
    # Fails when given a value that isn't a character vector
    expect_error(
        zscore_table_to_corecmotifs(
            example_zscore_table,
            zscore_columns = 8
        ),
        "zscore_columns is not a character vector"
    )

    # Fails when given a character vector that doesn't match the column names
    expect_error(
        zscore_table_to_corecmotifs(
            example_zscore_table,
            zscore_columns = c("cond1", "cond2")
        ),
        "zscore_table is missing one or more expected columns"
    )
})

test_that("bad arguments for probe_set_name are handled correctly", {
    # Fails when given a value that isn't a character vector
    expect_error(
        make_zscore_motif(
            example_zscore_table,
            probe_set_name = NA,
            pbm_condition = "UT_SUDHL4_SWISNF_mix"
        ),
        "probe_set_name is not a string"
    )

    # Fails when given a character vector with length > 1
    expect_error(
        make_zscore_motif(
            example_zscore_table,
            probe_set_name = c("MA0079.3_SP1", "another_probe_set"),
            pbm_condition = "UT_SUDHL4_SWISNF_mix"
        ),
        "probe_set_name is not a string"
    )
})

test_that("bad arguments for pbm_condition are handled correctly", {
    # Fails when given a value that isn't a character vector
    expect_error(
        make_zscore_motif(
            example_zscore_table,
            probe_set_name = "MA0079.3_SP1",
            pbm_condition = TRUE
        ),
        "pbm_condition is not a string"
    )

    # Fails when given a character vector with length > 1
    expect_error(
        make_zscore_motif(
            example_zscore_table,
            probe_set_name = "MA0079.3_SP1",
            pbm_condition = c("UT_SUDHL4_SWISNF_mix", ":o")
        ),
        "pbm_condition is not a string"
    )
})

test_that("bad arguments for array_id are handled correctly", {
    # Fails when given a value that isn't a character vector
    expect_error(
        zscore_table_to_corecmotifs(
            example_zscore_table,
            zscore_columns = colnames(example_fluorescence_table)[2:5],
            array_id = 4
        ),
        "array_id is not a string"
    )

    # Fails when given a character vector with length > 1
    expect_error(
        zscore_table_to_corecmotifs(
            example_zscore_table,
            zscore_columns = colnames(example_fluorescence_table)[2:5],
            array_id = c("too", "many", "strings")
        ),
        "array_id is not a string"
    )
})

test_that("example CoRecMotifs are reproducible", {
    # Convert the example z-score table into CoRecMotifs
    corecmotifs_1 <-
        zscore_table_to_corecmotifs(
            example_zscore_table,
            zscore_columns = colnames(example_fluorescence_table)[2:5],
            array_id = "v1_a11_run1"
        )

    # Make sure they're equal to the example CoRecMotifs
    expect_equal(corecmotifs_1, example_corecmotifs[1:20])

    # Convert one probe set/PBM condition combo into a z-score motif
    zscore_motif_1 <-
        make_zscore_motif(
            example_zscore_table,
            probe_set_name = "MA0079.3_SP1",
            pbm_condition = "UT_SUDHL4_SWISNF_mix"
        )

    # Make sure it's equal to the example CoRecMotif's z-score motif
    expect_equal(zscore_motif_1, example_corecmotifs[[1]]@zscore_motif)
})

