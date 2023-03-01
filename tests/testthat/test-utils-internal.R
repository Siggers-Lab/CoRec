test_that("update_output_base_name works()", {
    test_dir <- system.file("extdata", package = "CoRec")
    expect_equal(
        update_output_base_name(),
        NULL
    )
    expect_equal(
        update_output_base_name(test_dir),
        paste0(test_dir, "/output")
    )
    expect_equal(
        update_output_base_name(test_dir, output_base_name = "base_name"),
        paste0(test_dir, "/base_name")
    )
    expect_equal(
        update_output_base_name(test_dir, array_id = "sup"),
        paste0(test_dir, "/output_sup")
    )
    expect_equal(
        update_output_base_name(
            test_dir, output_base_name = "hello", array_id = "world"
        ),
        paste0(test_dir, "/hello_world")
    )
    expect_equal(
        update_output_base_name(output_base_name = "base_name"),
        paste0(getwd(), "/base_name")
    )
    expect_equal(
        update_output_base_name(output_base_name = "hello", array_id = "world"),
        paste0(getwd(), "/hello_world")
    )
    expect_equal(
        update_output_base_name(array_id = "bye"),
        NULL
    )
})

test_that("check_colnames works()", {
    # Make a test data frame
    test_df <-
        data.frame(
            "col_1" = rep(c(1, 2), 5),
            "col_2" = c(rep("a", 5), rep("b", 5)),
            "col_3" = letters[1:10]
        )

    # Fails when expected column names are missing
    expect_error(
        check_colnames(test_df, c("col_1", "col_2", "col_3", "col_4")),
        "test_df is missing one or more expected columns"
    )

    # Extra columns are removed
    expect_equal(
        check_colnames(test_df, c("col_2", "col_3")),
        data.frame(
            "col_2" = c(rep("a", 5), rep("b", 5)),
            "col_3" = letters[1:10]
        )
    )

    # Duplicate rows are removed
    expect_equal(
        check_colnames(test_df, c("col_1", "col_2")),
        data.frame(
            "col_1" = c(1, 2, 2, 1),
            "col_2" = c(rep("a", 2), rep("b", 2))
        )
    )
})

test_that("check_valid_zscore_motif() works", {
    # Fails when given a value that is numeric but isn't a matrix
    expect_error(
        check_valid_zscore_motif(c(0.1, 0.2, 0.3, 0.4)),
        "zscore_motif is not a matrix"
    )

    # Fails when given a value that is a matrix but isn't numeric
    expect_error(
        check_valid_zscore_motif(matrix(c("a", "b", "c", "d"))),
        "zscore_motif is not a numeric"
    )

    # Fails when given a numeric matrix with the wrong dimensions
    expect_error(
        check_valid_zscore_motif(matrix(rep(1, 20), nrow = 2)),
        "zscore_motif must have four rows"
    )

    # Gives a warning when given a numeric matrix with the wrong row names
    expect_warning(
        check_valid_zscore_motif(matrix(rep(1, 20), nrow = 4)),
        "zscore_motif row names are not 'A', 'C', 'G', 'T'"
    )

    # Make an example z-score motif with row names
    zscore_motif <- matrix(rep(1, 20), nrow = 4)
    rownames(zscore_motif) <- c("A", "C", "G", "T")

    expected_zscore_motif <- zscore_motif
    colnames(expected_zscore_motif) <-
        as.character(1:ncol(expected_zscore_motif))

    expect_equal(
        check_valid_zscore_motif(zscore_motif), expected_zscore_motif
    )
})

test_that("try_catch_save_output() works", {
    # Fails when given a value that isn't a character vector
    expect_error(
        try_catch_save_output(
            example_corecmotif_list_1,
            output_file = 4
        ),
        "output_file is not a string"
    )

    # Fails when given a character vector with length > 1
    expect_error(
        try_catch_save_output(
            example_corecmotif_list_1,
            output_file = c("too", "many", "strings")
        ),
        "output_file is not a string"
    )

    # Gives a warning when given a path to a file that can't be created
    expect_warning(
        try_catch_save_output(
            example_corecmotif_list_1,
            output_file = paste0(
                system.file(
                    "data/",
                    package = "CoRec"
                ),
                "/nonexistent_folder/nonexistent_file.txt"
            )
        ),
        "Could not write to output file .*Skipping output file creation..."
    )
})

test_that("create_array_id() works", {
    withr::with_seed(
        seed = 4023488,
        expect_equal(create_array_id(), "random_id_73486761")
    )
})

test_that("find_seed_zscore() works", {
    zscore_motif <-
        matrix(
            c(
                c(0.5, 2.0, 0.5, 0.2, 0.2, 1.3),
                c(1.0, 0.5, 1.1, 0.3, 1.3, 0.5),
                c(1.2, 0.3, 1.2, 0.4, 1.5, 0.2),
                c(0.3, 0.4, 0.1, 0.5, 0.5, 0.3)
            ),
            nrow = 4,
            byrow = TRUE
        )

    expect_equal(find_seed_zscore(zscore_motif), 0.5)
})

test_that("calculate_beta() works", {
    zscore_motif_1 <-
        matrix(
            c(
                c(0.5, 2.0, 0.5, 0.2, 0.2, 1.3),
                c(1.0, 0.5, 1.1, 0.3, 1.3, 0.5),
                c(1.2, 0.3, 1.2, 0.4, 1.5, 0.2),
                c(0.3, 0.4, 0.1, 0.5, 0.5, 0.3)
            ),
            nrow = 4,
            byrow = TRUE
        )

    expect_equal(calculate_beta(zscore_motif_1), 3.75)

    zscore_motif_2 <-
        matrix(
            c(
                c(-0.5, 2.0, -0.5, 0.2, 0.2, 1.3),
                c(1.0, -0.5, 1.1, 0.3, 1.3, -0.5),
                c(1.2, 0.3, 1.2, 0.4, 1.5, 0.2),
                c(0.3, 0.4, 0.1, -0.5, -0.5, 0.3)
            ),
            nrow = 4,
            byrow = TRUE
        )

    expect_equal(calculate_beta(zscore_motif_2), 4)

    zscore_motif_3 <-
        matrix(
            c(
                c(8.0, 2.0, 8.0, 0.2, 0.2, 1.3),
                c(1.0, 8.0, 1.1, 0.3, 1.3, 8.0),
                c(1.2, 0.3, 1.2, 0.4, 1.5, 0.2),
                c(0.3, 0.4, 0.1, 8.0, 8.0, 0.3)
            ),
            nrow = 4,
            byrow = TRUE
        )

    expect_equal(calculate_beta(zscore_motif_3), 1)
})

test_that("calculate_strength() works", {
    zscore_motif <-
        matrix(
            c(
                c(0.50, 2.02, 0.50, 0.24, 0.25, 1.36),
                c(1.01, 0.50, 1.13, 0.34, 1.35, 0.50),
                c(1.21, 0.32, 1.23, 0.44, 1.55, 0.26),
                c(0.31, 0.42, 0.13, 0.50, 0.50, 0.36)
            ),
            nrow = 4,
            byrow = TRUE
        )

    expect_equal(calculate_strength(zscore_motif), 1.5006)
})

test_that("zscore_to_universalmotif() works", {
    expect_equal(
        zscore_to_universalmotif(
            get_zscore_motif(example_corecmotifs[[6]]),
            beta = get_beta(example_corecmotifs[[6]]),
            name = get_motif_name(example_corecmotifs[[6]])
        ),
        get_motif(example_corecmotifs[[6]])
    )
})

test_that("calculate_rolling_ic() works", {
    expect_equal(
        calculate_rolling_ic(get_motif(example_corecmotifs[[15]])),
        get_rolling_ic(example_corecmotifs[[15]])
    )
})

test_that("sum_negatives() and sum_positives() work", {
    motif_matrix <-
        matrix(
            c(
                c(-1.0, 2.00, 3.00, 4.00, 1.00),
                c(-2.0, -3.0, 2.00, 3.00, 1.00),
                c(-1.0, -4.0, -3.0, 2.00, 3.00),
                c(-4.0, -5.0, -4.0, -2.0, 5.00)
            ),
            nrow = 4,
            byrow = TRUE
        )

    expect_equal(sum_negatives(motif_matrix), c(-8.0, -12.0, -7.0, -2.0, 0))
    expect_equal(sum_positives(motif_matrix), c(0.0, 2.0, 5.0, 9.0, 10.0))
})

