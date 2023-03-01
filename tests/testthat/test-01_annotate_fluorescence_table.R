test_that("bad arguments for fluorescence_table are handled correctly", {
    # Fails when given a value that isn't a data frame
    expect_error(
        annotate_fluorescence_table(
            4,
            fluorescence_columns = c(
                "UT_SUDHL4_SWISNF_mix", "UT_SUDHL4_HDAC_mix"
            ),
            annotation = hTF_v1_annotation
        ),
        "fluorescence_table is not a data frame"
    )

    # Fails when given a data frame with the wrong column names
    expect_error(
        annotate_fluorescence_table(
            data.frame("probe_id" = c("1", "2", "3"), "column_2" = c(4, 5, 6)),
            fluorescence_columns = c(
                "UT_SUDHL4_SWISNF_mix", "UT_SUDHL4_HDAC_mix"
            ),
            annotation = hTF_v1_annotation
        ),
        "fluorescence_table is missing one or more expected columns"
    )
})

test_that("bad arguments for fluorescence_columns are handled correctly", {
    # Fails when given a value that isn't a character vector
    expect_error(
        annotate_fluorescence_table(
            example_fluorescence_table,
            fluorescence_columns = 8,
            annotation = hTF_v1_annotation
        ),
        "fluorescence_columns is not a character vector"
    )

    # Fails when given a character vector that doesn't match the column names
    expect_error(
        annotate_fluorescence_table(
            example_fluorescence_table,
            fluorescence_columns = c("cond1", "cond2"),
            annotation = hTF_v1_annotation
        ),
        "fluorescence_table is missing one or more expected columns"
    )
})

test_that("bad arguments for annotation are handled correctly", {
    # Fails when given a value that isn't a data frame
    expect_error(
        annotate_fluorescence_table(
            example_fluorescence_table,
            fluorescence_columns = colnames(example_fluorescence_table)[2:5],
            annotation = c("this", "is", "not", "a", "data", "frame")
        ),
        "annotation is not a data frame"
    )

    # Fails when given a data frame with the wrong column names
    expect_error(
        annotate_fluorescence_table(
            example_fluorescence_table,
            fluorescence_columns = colnames(example_fluorescence_table)[2:5],
            annotation = data.frame(
                "probe_id" = c("1", "2", "3"),
                "column_2" = c(4, 5, 6)
            ),
        ),
        "annotation is missing one or more expected columns"
    )

    # Gives a warning when annotation is missing probe IDs in fluorescence_table
    expect_warning(
        annotate_fluorescence_table(
            example_fluorescence_table,
            fluorescence_columns = colnames(example_fluorescence_table)[2:5],
            annotation = dplyr::filter(
                hTF_v1_annotation,
                probe_set %in% c("MA0785.1_POU2F1", "MA0686.1_SPDEF")
            ),
        ),
        "annotation is missing probe IDs present in fluorescence_table"
    )
})

test_that("example annotated fluorescence tables are reproducible", {
    # Annotate the example un-annotated fluorescence table
    annotated_fluorescence_1 <-
        annotate_fluorescence_table(
            example_fluorescence_table,
            fluorescence_columns = colnames(example_fluorescence_table)[2:5],
            annotation = hTF_v1_annotation
        )

    # Make sure it's equal to the example annotated fluorescence table
    expect_equal(annotated_fluorescence_1, example_annotated_fluorescence_table)

    # Annotate part of the example un-annotated fluorescence table
    annotated_fluorescence_2 <-
        annotate_fluorescence_table(
            example_fluorescence_table,
            fluorescence_columns = colnames(example_fluorescence_table[2:3]),
            annotation = hTF_v1_annotation
        )

    # Make sure the extra columns got dropped
    expect_equal(
        annotated_fluorescence_2, example_annotated_fluorescence_table[1:8]
    )
})

