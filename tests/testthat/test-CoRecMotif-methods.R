test_that("as.data.frame() works", {
    # Make a data frame to compare against for an unmatched motif
    expected_df_1 <-
        data.frame(
            "motif_name" = example_corecmotifs[[1]]@motif@name,
            "probe_set" = example_corecmotifs[[1]]@probe_set,
            "pbm_condition" = example_corecmotifs[[1]]@pbm_condition,
            "array_id" = example_corecmotifs[[1]]@array_id,
            "motif_strength" = example_corecmotifs[[1]]@motif_strength,
            "rolling_ic" = example_corecmotifs[[1]]@rolling_ic,
            "seed_sequence" = example_corecmotifs[[1]]@seed_sequence,
            "match_motif" = NA,
            "match_pvalue" = NA_real_,
            "match_cluster" = NA_character_
        )

    expect_equal(as.data.frame(example_corecmotifs[[1]]), expected_df_1)

    # Make a data frame to compare against
    expected_df_2 <-
        data.frame(
            "motif_name" = example_matched_corecmotifs[[7]]@motif@name,
            "probe_set" = example_matched_corecmotifs[[7]]@probe_set,
            "pbm_condition" = example_matched_corecmotifs[[7]]@pbm_condition,
            "array_id" = example_matched_corecmotifs[[7]]@array_id,
            "motif_strength" = example_matched_corecmotifs[[7]]@motif_strength,
            "rolling_ic" = example_matched_corecmotifs[[7]]@rolling_ic,
            "seed_sequence" = example_matched_corecmotifs[[7]]@seed_sequence,
            "match_motif" =
                example_matched_corecmotifs[[7]]@match_motif@altname,
            "match_pvalue" = example_matched_corecmotifs[[7]]@match_pvalue,
            "match_cluster" = example_matched_corecmotifs[[7]]@match_cluster
        )

    expect_equal(as.data.frame(example_matched_corecmotifs[[7]]), expected_df_2)
})

test_that("getters work", {
    # Get an example motif and calculate some expected values
    motif_1 <- example_corecmotifs[[1]]

    delta_zscore_1 <-
        motif_1@zscore_motif %>%
        as.data.frame() %>%
        dplyr::mutate_all(list(~ . - median(.))) %>%
        as.matrix()

    seed_zscore_1 <-
        Reduce(dplyr::intersect, as.list(as.data.frame(motif_1@zscore_motif)))

    icm_1 <- universalmotif::convert_type(motif_1@motif, "ICM")["motif"]
    pwm_1 <-
        suppressMessages(
            universalmotif::convert_type(motif_1@motif, "PWM")["motif"]
        )
    ppm_1 <- universalmotif::convert_type(motif_1@motif, "PPM")["motif"]

    # Make sure everything matches
    expect_equal(get_probe_set(motif_1), motif_1@probe_set)
    expect_equal(get_pbm_condition(motif_1), motif_1@pbm_condition)
    expect_equal(get_zscore_motif(motif_1), motif_1@zscore_motif)
    expect_equal(get_delta_zscore_motif(motif_1), delta_zscore_1)

    expect_equal(get_array_id(motif_1), motif_1@array_id)
    expect_equal(get_motif_strength(motif_1), motif_1@motif_strength)
    expect_equal(get_rolling_ic(motif_1), motif_1@rolling_ic)
    expect_equal(get_seed_zscore(motif_1), seed_zscore_1)
    expect_equal(get_seed_sequence(motif_1), motif_1@seed_sequence)
    expect_equal(get_beta(motif_1), motif_1@beta)

    expect_equal(get_motif(motif_1), motif_1@motif)
    expect_equal(get_motif_name(motif_1), motif_1@motif@name)
    expect_equal(get_icm(motif_1), icm_1)
    expect_equal(suppressMessages(get_pwm(motif_1)), pwm_1)
    expect_equal(get_ppm(motif_1), ppm_1)

    expect_equal(get_match_motif(motif_1), motif_1@match_motif)
    expect_equal(get_match_name(motif_1), NA)
    expect_equal(get_match_altname(motif_1), NA)
    expect_equal(get_match_icm(motif_1), NA)
    expect_equal(get_match_pwm(motif_1), NA)
    expect_equal(get_match_ppm(motif_1), NA)

    expect_equal(get_match_pvalue(motif_1), motif_1@match_pvalue)
    expect_equal(get_match_cluster(motif_1), motif_1@match_cluster)

    # Get an example motif with a match and calculate some expected values
    motif_2 <- example_matched_corecmotifs[[3]]

    icm_2 <- universalmotif::convert_type(motif_2@match_motif, "ICM")
    icm_2_rc <- universalmotif::motif_rc(icm_2)
    pwm_2 <- universalmotif::convert_type(motif_2@match_motif, "PWM")
    pwm_2_rc <- universalmotif::motif_rc(pwm_2)
    ppm_2 <- universalmotif::convert_type(motif_2@match_motif, "PPM")
    ppm_2_rc <- universalmotif::motif_rc(ppm_2)

    # Make sure everything matches
    expect_equal(get_match_motif(motif_2), motif_2@match_motif)
    expect_equal(get_match_name(motif_2), motif_2@match_motif["name"])
    expect_equal(get_match_altname(motif_2), motif_2@match_motif["altname"])
    expect_equal(
        get_match_icm(motif_2, correct_orientation = FALSE), icm_2["motif"]
    )
    expect_equal(
        get_match_pwm(motif_2, correct_orientation = FALSE), pwm_2["motif"]
    )
    expect_equal(
        get_match_ppm(motif_2, correct_orientation = FALSE), ppm_2["motif"]
    )
    expect_equal(
        get_match_icm(motif_2, correct_orientation = TRUE), icm_2_rc["motif"]
    )
    expect_equal(
        get_match_pwm(motif_2, correct_orientation = TRUE), pwm_2_rc["motif"]
    )
    expect_equal(
        get_match_ppm(motif_2, correct_orientation = TRUE), ppm_2_rc["motif"]
    )
})

test_that("bad arguments for character field setters are handled correctly", {
    # All fail when given a value that isn't a character vector
    expect_error(
        set_probe_set(example_corecmotifs[[11]], 4),
        "assignment of an object of class .* is not valid"
    )

    expect_error(
        set_pbm_condition(example_corecmotifs[[18]], NA),
        "assignment of an object of class .* is not valid"
    )

    expect_error(
        set_array_id(example_corecmotifs[[35]], NULL),
        "assignment of an object of class .* is not valid"
    )

    expect_error(
        set_seed_sequence(example_corecmotifs[[6]], data.frame()),
        "assignment of an object of class .* is not valid"
    )

    expect_error(
        set_motif_name(example_corecmotifs[[8]], list()),
        "assignment of an object of class .* is not valid"
    )

    # All fail when given a character vector with length > 1
    expect_error(
        set_probe_set(example_corecmotifs[[10]], c("aaaaaa", "help")),
        "@probe_set must be a character vector of length 1"
    )

    expect_error(
        set_pbm_condition(example_corecmotifs[[12]], c("aaaaaa", "help")),
        "@pbm_condition must be a character vector of length 1"
    )

    expect_error(
        set_array_id(example_corecmotifs[[26]], c("I'm", "so", "tired")),
        "@array_id must be a character vector of length 1"
    )

    expect_error(
        set_seed_sequence(example_corecmotifs[[3]], c("aaaaaa", "zzzzzz")),
        "@seed_sequence must be a character vector of length 1"
    )

    expect_error(
        set_motif_name(example_corecmotifs[[19]], c("aaaaaa", ":'(")),
        "name must be length 1"
    )
})

test_that("character field setters work", {
    # Update @probe_set
    motif_1 <- example_corecmotifs[[21]]
    motif_1@probe_set <- "something new"

    expect_equal(
        set_probe_set(example_corecmotifs[[21]], "something new"),
        motif_1
    )

    # Update @pbm_condition
    motif_2 <- example_corecmotifs[[40]]
    motif_2@pbm_condition <- "blah"

    expect_equal(
        set_pbm_condition(example_corecmotifs[[40]], "blah"),
        motif_2
    )

    # Update @array_id
    motif_3 <- example_corecmotifs[[33]]
    motif_3@array_id <- "hello_world"

    expect_equal(
        set_array_id(example_corecmotifs[[33]], "hello_world"),
        motif_3
    )

    # Update @seed_sequence
    motif_4 <- example_corecmotifs[[9]]
    motif_4@array_id <- "AAAAAAAAAAAAAAA"

    expect_equal(
        set_array_id(example_corecmotifs[[9]], "AAAAAAAAAAAAAAA"),
        motif_4
    )

    # Update @motif@name
    motif_5 <- example_corecmotifs[[13]]
    motif_5@motif@name <- "bob"

    expect_equal(
        set_motif_name(example_corecmotifs[[13]], "bob"),
        motif_5
    )
})

test_that("bad arguments for set_zscore_motif() are handled correctly", {
    # Fails when given a value that is numeric but isn't a matrix
    expect_error(
        set_zscore_motif(
            example_corecmotifs[[22]],
            c(0.1, 0.2, 0.3, 0.4)
        ),
        "zscore_motif is not a matrix"
    )

    # Fails when given a value that is a matrix but isn't numeric
    expect_error(
        set_zscore_motif(
            example_corecmotifs[[23]],
            matrix(c("a", "b", "c", "d"))
        ),
        "zscore_motif is not a numeric"
    )

    # Fails when given a numeric matrix with the wrong dimensions
    expect_error(
        set_zscore_motif(
            example_corecmotifs[[24]],
            matrix(rep(1, 20), nrow = 2)
        ),
        "zscore_motif must have four rows"
    )

    # Gives a warning when given a numeric matrix with the wrong row names
    expect_warning(
        set_zscore_motif(
            example_corecmotifs[[25]],
            matrix(rep(1, 20), nrow = 4)
        ),
        "zscore_motif row names are not 'A', 'C', 'G', 'T'"
    )
})

test_that("set_zscore_motif() works", {
    # Update the z-score motif of an example CoRecMotif
    zscore_motif_1 <- get_zscore_motif(example_corecmotifs[[10]])
    motif_1 <- set_zscore_motif(example_corecmotifs[[20]], zscore_motif_1)

    # Make sure all the dependent slots got updated
    expect_equal(
        get_zscore_motif(motif_1),
        get_zscore_motif(example_corecmotifs[[10]])
    )
    expect_equal(
        get_motif_strength(motif_1),
        get_motif_strength(example_corecmotifs[[10]])
    )
    expect_equal(
        get_beta(motif_1),
        get_beta(example_corecmotifs[[10]])
    )
    expect_equal(
        get_ppm(motif_1),
        get_ppm(example_corecmotifs[[10]])
    )
    expect_equal(
        get_rolling_ic(motif_1),
        get_rolling_ic(example_corecmotifs[[10]])
    )

    # Update the z-score motif of an example CoRecMotif with a match
    motif_2 <-
        set_zscore_motif(example_matched_corecmotifs[[3]], zscore_motif_1)

    # Make sure all the match slots got updated
    expect_equal(get_match_motif(motif_2), NA)
    expect_equal(get_match_name(motif_2), NA)
    expect_equal(get_match_altname(motif_2), NA)
    expect_equal(get_match_icm(motif_2), NA)
    expect_equal(get_match_pwm(motif_2), NA)
    expect_equal(get_match_ppm(motif_2), NA)
    expect_equal(get_match_pvalue(motif_2), NA_real_)
    expect_equal(get_match_cluster(motif_2), NA_character_)
})

