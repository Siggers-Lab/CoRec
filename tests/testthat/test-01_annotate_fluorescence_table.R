test_that("bad arguments for fluorescence_file are handled correctly", {
    # Define some valid values for the other parameters
    test_pbm_conditions <- c("cond1", "cond2", "cond3", "cond4")
    test_annotation_file <- system.file(
        "example_data/hTF_v1_example_annotation.tsv",
        package = "hTFArrayAnalysis"
    )

    # Fails when given a value that isn't a character vector
    expect_error(
        annotate_fluorescence_table(
            4,
            pbm_conditions = test_pbm_conditions,
            annotation_file = test_annotation_file
        ),
        "fluorescence_file is not a string"
    )

    # Fails when given a character vector with length > 1
    expect_error(
        annotate_fluorescence_table(
            c(
                system.file(
                    "example_data/hTF_v1_example_fluorescence_rep1.dat",
                    package = "hTFArrayAnalysis"
                ),
                system.file(
                    "example_data/hTF_v1_example_fluorescence_rep2.dat",
                    package = "hTFArrayAnalysis"
                )
            ),
            pbm_conditions = test_pbm_conditions,
            annotation_file = test_annotation_file
        ),
        "fluorescence_file is not a string"
    )

    # Fails when given a path to a file that doesn't exist
    expect_error(
        annotate_fluorescence_table(
            system.file(
                "example_data/nonexistent_file.txt",
                package = "hTFArrayAnalysis"
            ),
            pbm_conditions = test_pbm_conditions,
            annotation_file = test_annotation_file
        ),
        "Path .* does not exist"
    )
})

test_that("bad arguments for pbm_conditions are handled correctly", {
    # Define some valid values for the other parameters
    test_fluorescence_file <- system.file(
        "example_data/hTF_v1_example_fluorescence_rep1.dat",
        package = "hTFArrayAnalysis"
    )
    test_annotation_file <- system.file(
        "example_data/hTF_v1_example_annotation.tsv",
        package = "hTFArrayAnalysis"
    )

    # Fails when given a value that isn't a character vector
    expect_error(
        annotate_fluorescence_table(
            test_fluorescence_file,
            pbm_conditions = 8,
            annotation_file = test_annotation_file
        ),
        "pbm_conditions is not a character vector"
    )

    # Fails when given a character vector with the wrong length
    expect_error(
        annotate_fluorescence_table(
            test_fluorescence_file,
            pbm_conditions = c("cond1", "cond2"),
            annotation_file = test_annotation_file
        ),
        "pbm_conditions is the wrong length"
    )
})

test_that("bad arguments for annotation_file are handled correctly", {
    # Define some valid values for the other parameters
    test_fluorescence_file <- system.file(
        "example_data/hTF_v1_example_fluorescence_rep1.dat",
        package = "hTFArrayAnalysis"
    )
    test_pbm_conditions <- c("cond1", "cond2", "cond3", "cond4")

    # Fails when given a value that isn't a character vector
    expect_error(
        annotate_fluorescence_table(
            test_fluorescence_file,
            pbm_conditions = test_pbm_conditions,
            annotation_file = 9
        ),
        "annotation_file is not a string"
    )

    # Fails when given a character vector with length > 1
    expect_error(
        annotate_fluorescence_table(
            test_fluorescence_file,
            pbm_conditions = test_pbm_conditions,
            annotation_file = c(
                system.file(
                    "example_data/hTF_v1_example_fluorescence_rep1.dat",
                    package = "hTFArrayAnalysis"
                ),
                system.file(
                    "example_data/hTF_v1_example_fluorescence_rep2.dat",
                    package = "hTFArrayAnalysis"
                )
            )
        ),
        "annotation_file is not a string"
    )

    # Fails when given a path to a file that doesn't exist
    expect_error(
        annotate_fluorescence_table(
            test_fluorescence_file,
            pbm_conditions = test_pbm_conditions,
            annotation_file = system.file(
                "example_data/nonexistent_file.txt",
                package = "hTFArrayAnalysis"
            )
        ),
        "Path .* does not exist"
    )

    # Fails when given a path to a file missing the expected columns
    expect_error(
        annotate_fluorescence_table(
            test_fluorescence_file,
            pbm_conditions = test_pbm_conditions,
            annotation_file = test_fluorescence_file
        ),
        "annotation_file is missing one or more expected columns"
    )
})

test_that("bad arguments for output_file are handled correctly", {
    # Define some valid values for the other parameters
    test_fluorescence_file <- system.file(
        "example_data/hTF_v1_example_fluorescence_rep1.dat",
        package = "hTFArrayAnalysis"
    )
    test_pbm_conditions <- c("cond1", "cond2", "cond3", "cond4")
    test_annotation_file <- system.file(
        "example_data/hTF_v1_example_annotation.tsv",
        package = "hTFArrayAnalysis"
    )

    # Fails when given a value that isn't a character vector
    expect_error(
        annotate_fluorescence_table(
            test_fluorescence_file,
            pbm_conditions = test_pbm_conditions,
            annotation_file = test_annotation_file,
            output_file = 4
        ),
        "output_file is not a string"
    )

    # Fails when given a character vector with length > 1
    expect_error(
        annotate_fluorescence_table(
            test_fluorescence_file,
            pbm_conditions = test_pbm_conditions,
            annotation_file = test_annotation_file,
            output_file = c(
                system.file(
                    "example_data/hTF_v1_example_fluorescence_rep1.dat",
                    package = "hTFArrayAnalysis"
                ),
                system.file(
                    "example_data/hTF_v1_example_fluorescence_rep2.dat",
                    package = "hTFArrayAnalysis"
                )
            )
        ),
        "output_file is not a string"
    )

    # Gives a warning when given a path to a file that can't be created
    expect_warning(
        annotate_fluorescence_table(
            test_fluorescence_file,
            pbm_conditions = test_pbm_conditions,
            annotation_file = test_annotation_file,
            output_file = paste0(
                system.file(
                    "example_data/output",
                    package = "hTFArrayAnalysis"
                ),
                "/nonexistent_folder/nonexistent_file.txt"
            )
        ),
        "Could not write to output file .*Skipping output file creation..."
    )
})

