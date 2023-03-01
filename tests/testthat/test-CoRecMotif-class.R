test_that("bad arguments for probe_set are handled correctly", {
    # Make a test z-score matrix
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
    rownames(zscore_motif) <- c("A", "C", "G", "T")

    # Fails when given a value that isn't a character vector
    expect_error(
        CoRecMotif(
            probe_set = NA,
            pbm_condition = "hello",
            zscore_motif = zscore_motif
        ),
        "probe_set is not a string"
    )

    # Fails when given a character vector with length > 1
    expect_error(
        CoRecMotif(
            probe_set = c("hello", "world"),
            pbm_condition = "hello",
            zscore_motif = zscore_motif
        ),
        "probe_set is not a string"
    )
})

test_that("bad arguments for pbm_condition are handled correctly", {
    # Make a test z-score matrix
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
    rownames(zscore_motif) <- c("A", "C", "G", "T")

    # Fails when given a value that isn't a character vector
    expect_error(
        CoRecMotif(
            probe_set = "hello",
            pbm_condition = 4,
            zscore_motif = zscore_motif
        ),
        "pbm_condition is not a string"
    )

    # Fails when given a character vector with length > 1
    expect_error(
        CoRecMotif(
            probe_set = "hello",
            pbm_condition = c("hi", "everyone"),
            zscore_motif = zscore_motif
        ),
        "pbm_condition is not a string"
    )
})

test_that("bad arguments for zscore_motif are handled correctly", {
    # Fails when given a value that is numeric but isn't a matrix
    expect_error(
        CoRecMotif(
            probe_set = "probe_set_1",
            pbm_condition = "pbm_condition_1",
            zscore_motif = c(0.1, 0.2, 0.3, 0.4)
        ),
        "zscore_motif is not a matrix"
    )

    # Fails when given a value that is a matrix but isn't numeric
    expect_error(
        CoRecMotif(
            probe_set = "probe_set_1",
            pbm_condition = "pbm_condition_1",
            zscore_motif = matrix(c("a", "b", "c", "d"))
        ),
        "zscore_motif is not a numeric"
    )

    # Fails when given a numeric matrix with the wrong dimensions
    expect_error(
        CoRecMotif(
            probe_set = "probe_set_1",
            pbm_condition = "pbm_condition_1",
            zscore_motif = matrix(rep(1, 20), nrow = 2)
        ),
        "zscore_motif must have four rows"
    )

    # Gives a warning when given a numeric matrix with the wrong row names
    expect_warning(
        CoRecMotif(
            probe_set = "probe_set_1",
            pbm_condition = "pbm_condition_1",
            zscore_motif = matrix(rep(1, 20), nrow = 4)
        ),
        "zscore_motif row names are not 'A', 'C', 'G', 'T'"
    )
})

test_that("bad arguments for array_id are handled correctly", {
    # Make a test z-score matrix
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
    rownames(zscore_motif) <- c("A", "C", "G", "T")

    # Fails when given a value that isn't a character vector
    expect_error(
        CoRecMotif(
            probe_set = "hello",
            pbm_condition = "hello",
            zscore_motif = zscore_motif,
            array_id = NA
        ),
        "array_id is not a string"
    )

    # Fails when given a character vector with length > 1
    expect_error(
        CoRecMotif(
            probe_set = "hello",
            pbm_condition = "hello",
            zscore_motif = zscore_motif,
            array_id = c("too", "many", "IDs")
        ),
        "array_id is not a string"
    )
})

test_that("bad arguments for seed_sequence are handled correctly", {
    # Make a test z-score matrix
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
    rownames(zscore_motif) <- c("A", "C", "G", "T")

    # Fails when given a value that isn't a character vector
    expect_error(
        CoRecMotif(
            probe_set = "hello",
            pbm_condition = "hello",
            zscore_motif = zscore_motif,
            seed_sequence = NA
        ),
        "seed_sequence is not a string"
    )

    # Fails when given a character vector with length > 1
    expect_error(
        CoRecMotif(
            probe_set = "hello",
            pbm_condition = "hello",
            zscore_motif = zscore_motif,
            seed_sequence = c("too", "many", "IDs")
        ),
        "seed_sequence is not a string"
    )
})

test_that("example CoRecMotifs are reproducible", {
    # Remake an example CoRecMotif
    test_corecmotif <-
        CoRecMotif(
            probe_set = get_probe_set(example_corecmotifs[[1]]),
            pbm_condition = get_pbm_condition(example_corecmotifs[[1]]),
            zscore_motif = get_zscore_motif(example_corecmotifs[[1]]),
            array_id = get_array_id(example_corecmotifs[[1]]),
            seed_sequence = get_seed_sequence(example_corecmotifs[[1]])
        )

    expect_equal(test_corecmotif, example_corecmotifs[[1]])
})

