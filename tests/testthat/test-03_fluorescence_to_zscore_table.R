test_that("bad arguments for fluorescence_table are handled correctly", {
    # Fails when given a value that isn't a data frame
    expect_error(
        fluorescence_to_zscore_table(
            4,
            fluorescence_columns = c(
                "UT_SUDHL4_SWISNF_mix", "UT_SUDHL4_HDAC_mix"
            )
        ),
        "fluorescence_table is not a data frame"
    )

    # Fails when given a data frame with the wrong column names
    expect_error(
        fluorescence_to_zscore_table(
            data.frame("probe_id" = c("1", "2", "3"), "column_2" = c(4, 5, 6)),
            fluorescence_columns = c(
                "UT_SUDHL4_SWISNF_mix", "UT_SUDHL4_HDAC_mix"
            )
        ),
        "fluorescence_table is missing one or more expected columns"
    )
})

test_that("bad arguments for fluorescence_columns are handled correctly", {
    # Fails when given a value that isn't a character vector
    expect_error(
        fluorescence_to_zscore_table(
            example_annotated_fluorescence_table,
            fluorescence_columns = 8
        ),
        "fluorescence_columns is not a character vector"
    )

    # Fails when given a character vector that doesn't match the column names
    expect_error(
        fluorescence_to_zscore_table(
            example_annotated_fluorescence_table,
            fluorescence_columns = c("cond1", "cond2")
        ),
        "fluorescence_table is missing one or more expected columns"
    )
})

test_that("example z-score tables are reproducible", {
    # Convert the example fluorescence table into z-scores
    zscores_1 <-
        fluorescence_to_zscore_table(
            example_annotated_fluorescence_table,
            fluorescence_columns = colnames(example_fluorescence_table)[2:5]
        )

    # Make sure it's equal to the example z-score table
    expect_equal(zscores_1, example_zscore_table)

    # Convert part of the example fluorescence table into z-scores
    zscores_2 <-
        fluorescence_to_zscore_table(
            example_annotated_fluorescence_table,
            fluorescence_columns = colnames(example_fluorescence_table[2:3])
        )

    # Make sure the extra columns got dropped
    expect_equal(zscores_2, example_zscore_table[1:8])
})

