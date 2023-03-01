test_that("bad arguments for n_replicates are handled correctly", {
    # Fails when given a value that isn't a numeric vector
    expect_error(
        check_replicates(
            example_corecmotifs,
            n_replicates = "hello"
        ),
        "n_replicates is not a count"
    )

    # Fails when given a numeric vector with length > 1
    expect_error(
        check_replicates(
            example_corecmotifs,
            n_replicates = c(1, 2)
        ),
        "n_replicates is not a count"
    )

    # Fails when given a non-integer
    expect_error(
        check_replicates(
            example_corecmotifs,
            n_replicates = 1.2
        ),
        "n_replicates is not a count"
    )
})

test_that("bad arguments for eucl_distance are handled correctly", {
    # Fails when given a value that isn't a numeric vector
    expect_error(
        check_replicates(
            example_corecmotifs,
            eucl_distance = "hello"
        ),
        "eucl_distance is not a number"
    )

    # All fail when given a numeric vector with length > 1
    expect_error(
        check_replicates(
            example_corecmotifs,
            eucl_distance = c(1, 2)
        ),
        "eucl_distance is not a number"
    )
})

test_that("filtering by n_replicates works", {
    # Regroup some of the example CoRecMotifs for testing purposes
    test_corecmotifs <-
        example_corecmotifs[1:12] %>%
        lapply(set_pbm_condition, "test")

    # Filter out groups with fewer than 2 motifs
    n_replicates_filter_1 <-
        check_replicates(
            test_corecmotifs,
            n_replicates = 2,
            eucl_distance = NULL
        )

    expect_setequal(n_replicates_filter_1, test_corecmotifs)

    # Filter out groups with fewer than 5 motifs
    n_replicates_filter_1 <-
        check_replicates(
            test_corecmotifs,
            n_replicates = 5,
            eucl_distance = NULL
        )

    expect_setequal(n_replicates_filter_1, list())

    # Filter out groups with fewer than 2 motifs after subsetting
    n_replicates_filter_2 <-
        check_replicates(
            test_corecmotifs[c(1:2, 4, 8, 11)],
            n_replicates = 2,
            eucl_distance = NULL
        )

    expect_setequal(n_replicates_filter_2, test_corecmotifs[c(1:2, 4)])

    # Filter out groups with fewer than 3 motifs after subsetting
    n_replicates_filter_2 <-
        check_replicates(
            test_corecmotifs[c(1:4, 6:8, 11)],
            n_replicates = 3,
            eucl_distance = NULL
        )

    expect_setequal(n_replicates_filter_2, test_corecmotifs[c(1:4, 6:8)])
})

test_that("filtering by eucl_distance works", {
    # Regroup some of the example CoRecMotifs for testing purposes
    test_corecmotifs <-
        example_corecmotifs[1:12] %>%
        lapply(set_pbm_condition, "test")

    # Filter out motifs with no Euclidean distance less than 1000
    eucl_distance_filter_1 <-
        check_replicates(
            test_corecmotifs,
            n_replicates = 1,
            eucl_distance = 1000
        )

    expect_setequal(eucl_distance_filter_1, test_corecmotifs)

    # Filter out motifs with no Euclidean distance less than 0.6
    eucl_distance_filter_2 <-
        check_replicates(
            test_corecmotifs,
            n_replicates = 1,
            eucl_distance = 0.6
        )

    expect_setequal(eucl_distance_filter_2, test_corecmotifs[c(1, 3, 5:12)])

    # Filter out motifs with no Euclidean distance less than 0.4
    eucl_distance_filter_3 <-
        check_replicates(
            test_corecmotifs,
            n_replicates = 1,
            eucl_distance = 0.4
        )

    expect_setequal(eucl_distance_filter_3, test_corecmotifs[9:10])

    # Filter out motifs with no Euclidean distance less than 0.1
    eucl_distance_filter_4 <-
        check_replicates(
            test_corecmotifs,
            n_replicates = 1,
            eucl_distance = 0.1
        )

    expect_setequal(eucl_distance_filter_4, list())
})

test_that("filtering by both n_replicates and eucl_distance works", {
    # Regroup some of the example CoRecMotifs for testing purposes
    test_corecmotifs <-
        example_corecmotifs[1:12] %>%
        lapply(set_pbm_condition, "test")

    # Filter out groups with < 4 motifs and motifs with no ED < 1000
    combined_filter_1 <-
        check_replicates(
            test_corecmotifs,
            n_replicates = 4,
            eucl_distance = 1000
        )

    expect_setequal(combined_filter_1, test_corecmotifs)

    # Filter out groups with < 3 motifs and motifs with no ED < 0.6
    combined_filter_2 <-
        check_replicates(
            test_corecmotifs,
            n_replicates = 3,
            eucl_distance = 0.6
        )

    expect_setequal(combined_filter_2, test_corecmotifs[5:12])

    # Filter out groups with < 2 motifs and motifs with no ED < 0.4
    combined_filter_3 <-
        check_replicates(
            test_corecmotifs,
            n_replicates = 2,
            eucl_distance = 0.4
        )

    expect_setequal(combined_filter_3, test_corecmotifs[9:10])

    # Filter out groups with < 3 motifs and motifs with no ED < 0.4
    combined_filter_4 <-
        check_replicates(
            test_corecmotifs,
            n_replicates = 3,
            eucl_distance = 0.4
        )

    expect_setequal(combined_filter_4, list())

    # Filter out subsetted groups with < 3 motifs and motifs with no ED < 0.6
    combined_filter_5 <-
        check_replicates(
            test_corecmotifs[c(1:4, 6:8, 11)],
            n_replicates = 2,
            eucl_distance = 0.6
        )

    expect_setequal(combined_filter_5, test_corecmotifs[c(1, 3, 7:8)])
})

