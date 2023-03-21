test_that("bad arguments for fluorescence_file are handled correctly", {
    # Fails when given a value that isn't a character vector
    expect_error(
        load_fluorescence_data(
            4,
            pbm_conditions = c("1", "2", "3", "4")
        ),
        "fluorescence_file is not a string"
    )

    # Fails when given a character vector with length > 1
    expect_error(
        load_fluorescence_data(
            c("two", "strings"),
            pbm_conditions = c("1", "2", "3", "4")
        ),
        "fluorescence_file is not a string"
    )
})

test_that("bad arguments for pbm_conditions are handled correctly", {
    # Fails when given a value that isn't a character vector
    expect_error(
        load_fluorescence_data(
            system.file(
                "extdata/example_fluorescence_data_1.dat",
                package = "CoRec"
            ),
            pbm_conditions = NA
        ),
        "pbm_conditions is not a character vector"
    )

    # Fails when given a character vector with the wrong length
    expect_error(
        load_fluorescence_data(
            system.file(
                "extdata/example_fluorescence_data_1.dat",
                package = "CoRec"
            ),
            pbm_conditions = c("this", "is", "too", "many", "strings")
        ),
        "pbm_conditions is the wrong length"
    )
})

test_that("example fluorescence table is reproducible", {
    pbm_conditions <-
        c(
            "UT_SUDHL4_SWISNF_mix",
            "UT_SUDHL4_HDAC_mix",
            "UT_SUDHL4_PRMT5",
            "UT_SUDHL4_JMJD2A"
        )

    fluorescence_table <-
        load_fluorescence_data(
            system.file(
                "extdata/example_fluorescence_data_1.dat",
                package = "CoRec"
            ),
            pbm_conditions = pbm_conditions
        )

    expect_equal(fluorescence_table, example_fluorescence_table)
})

