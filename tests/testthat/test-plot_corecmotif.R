test_that("plot_corecmotif() works", {
    # Make sure all the CoRecMotif logo types work
    expect_s3_class(
        plot_corecmotif(
            example_matched_corecmotifs[[1]],
            corecmotif_logo_type = "delta_zscore"
        ),
        "gg"
    )
    expect_s3_class(
        plot_corecmotif(
            example_matched_corecmotifs[[1]],
            corecmotif_logo_type = "ICM"
        ),
        "gg"
    )
    expect_s3_class(
        suppressMessages(
            plot_corecmotif(
                example_matched_corecmotifs[[1]],
                corecmotif_logo_type = "PWM"
            )
        ),
        "gg"
    )
    expect_s3_class(
        plot_corecmotif(
            example_matched_corecmotifs[[1]],
            corecmotif_logo_type = "PPM"
        ),
        "gg"
    )
    expect_s3_class(
        plot_corecmotif(
            example_matched_corecmotifs[[1]],
            corecmotif_logo_type = "none"
        ),
        "gg"
    )

    # Make sure all the reference logo types work
    expect_s3_class(
        plot_corecmotif(
            example_matched_corecmotifs[[1]],
            reference_logo_type = "ICM"
        ),
        "gg"
    )
    expect_s3_class(
        plot_corecmotif(
            example_matched_corecmotifs[[1]],
            reference_logo_type = "PWM"
        ),
        "gg"
    )
    expect_s3_class(
        plot_corecmotif(
            example_matched_corecmotifs[[1]],
            reference_logo_type = "PPM"
        ),
        "gg"
    )
    expect_s3_class(
        plot_corecmotif(
            example_matched_corecmotifs[[1]],
            reference_logo_type = "none"
        ),
        "gg"
    )

    # Make sure nothing is returned if both plot types are set to none
    expect_equal(
        plot_corecmotif(
            example_matched_corecmotifs[[1]],
            corecmotif_logo_type = "none",
            reference_logo_type = "none"
        ),
        NULL
    )
    expect_equal(
        plot_corecmotif(
            example_corecmotifs[[1]],
            corecmotif_logo_type = "none",
            reference_logo_type = "ICM"
        ),
        NULL
    )
})
