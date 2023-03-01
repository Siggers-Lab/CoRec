test_that("bad arguments for reference_motifs_file are handled correctly", {
    skip_if(!memes::meme_is_installed(), "MEME is not installed")

    # Fails when given a value that isn't a character vector
    expect_error(
        find_match(
            example_corecmotifs[1:4],
            reference_motifs_file = 4
        ),
        "reference_motifs_file is not a string"
    )

    # Fails when given a character vector with length > 1
    expect_error(
        find_match(
            example_corecmotifs[1:4],
            reference_motifs_file = c("two", "strings")
        ),
        "reference_motifs_file is not a string"
    )

    # Fails when given a path that isn't a MEME file
    suppressMessages(expect_error(
        find_match(
            example_corecmotifs[1:4],
            reference_motifs_file =
                system.file(
                    "extdata/example_fluorescence_data_1.dat",
                    package = "CoRec"
                )
        ),
        "Shell process had non-zero exit status."
    ))
})

test_that("bad arguments for cluster_assignments are handled correctly", {
    skip_if(!memes::meme_is_installed(), "MEME is not installed")

    # Fails when given a value that isn't a data frame
    expect_error(
        find_match(
            example_corecmotifs[1:4],
            reference_motifs_file =
                system.file(
                    "extdata/Homo_sapiens_JASPAR2022_CORE_filtered.meme",
                    package = "CoRec"
                ),
            cluster_assignments =
                c("motif_1" = "cluster_1", "motif_2" = "cluster_2")
        ),
        "cluster_assignments is not a data frame"
    )

    # Fails when given a data frame with the wrong column names
    expect_error(
        find_match(
            example_corecmotifs[1:4],
            reference_motifs_file =
                system.file(
                    "extdata/Homo_sapiens_JASPAR2022_CORE_filtered.meme",
                    package = "CoRec"
                ),
            cluster_assignments =
                data.frame(
                    "V1" = c("motif_1", "motif_2"),
                    "V2" = c("cluster_1", "cluster_2")
                )
        ),
        "cluster_assignments is missing one or more expected columns"
    )
})

test_that("bad arguments for meme_path are handled correctly", {
    skip_if(!memes::meme_is_installed(), "MEME is not installed")

    # Fails when given a value that isn't a character vector
    expect_error(
        find_match(
            example_corecmotifs[1:4],
            reference_motifs_file =
                system.file(
                    "extdata/Homo_sapiens_JASPAR2022_CORE_filtered.meme",
                    package = "CoRec"
                ),
            meme_path = 4
        ),
        "meme_path is not a string"
    )

    # Fails when given a character vector with length > 1
    expect_error(
        find_match(
            example_corecmotifs[1:4],
            reference_motifs_file =
                system.file(
                    "extdata/Homo_sapiens_JASPAR2022_CORE_filtered.meme",
                    package = "CoRec"
                ),
            meme_path = c("two", "strings")
        ),
        "meme_path is not a string"
    )

    # Fails when given a path that isn't to the MEME installation
    suppressMessages(expect_error(
        find_match(
            example_corecmotifs[1:4],
            reference_motifs_file =
                system.file(
                    "extdata/Homo_sapiens_JASPAR2022_CORE_filtered.meme",
                    package = "CoRec"
                ),
            meme_path =
                system.file(
                    "extdata/example_fluorescence_data_1.dat",
                    package = "CoRec"
                )
        ),
        "Invalid meme_path: could not find MEME installation."
    ))
})

test_that("bad arguments for min_overlap are handled correctly", {
    skip_if(!memes::meme_is_installed(), "MEME is not installed")

    # Fails when given a value that isn't a numeric vector
    expect_error(
        find_match(
            example_corecmotifs[1:4],
            reference_motifs_file =
                system.file(
                    "extdata/Homo_sapiens_JASPAR2022_CORE_filtered.meme",
                    package = "CoRec"
                ),
            min_overlap = "motif_1"
        ),
        "min_overlap is not a count"
    )

    # Fails when given a numeric vector with length > 1
    expect_error(
        find_match(
            example_corecmotifs[1:4],
            reference_motifs_file =
                system.file(
                    "extdata/Homo_sapiens_JASPAR2022_CORE_filtered.meme",
                    package = "CoRec"
                ),
            min_overlap = c(1, 2, 3)
        ),
        "min_overlap is not a count"
    )

    # Fails when given a non-integer
    expect_error(
        find_match(
            example_corecmotifs[1:4],
            reference_motifs_file =
                system.file(
                    "extdata/Homo_sapiens_JASPAR2022_CORE_filtered.meme",
                    package = "CoRec"
                ),
            min_overlap = 0.4
        ),
        "min_overlap is not a count"
    )
})

test_that("matching works", {
    skip_if(!memes::meme_is_installed(), "MEME is not installed")

    # Compare a CoRecMotif to the reference motifs and assign a cluster
    matched_motif_1 <-
        find_match(
            example_corecmotifs[8],
            reference_motifs_file =
                system.file(
                    "extdata/Homo_sapiens_JASPAR2022_CORE_filtered.meme",
                    package = "CoRec"
                ),
            cluster_assignments = motif_clusters
        )

    expect_equal(matched_motif_1, example_matched_corecmotifs[1])

    # Compare a CoRecMotif to the reference motifs without assigning a cluster
    matched_motif_2 <-
        find_match(
            example_corecmotifs[14],
            reference_motifs_file =
                system.file(
                    "extdata/Homo_sapiens_JASPAR2022_CORE_filtered.meme",
                    package = "CoRec"
                ),
        )

    # Remove the cluster assignment from the corresponding example motif
    matched_motif_no_cluster <- example_matched_corecmotifs[5]
    matched_motif_no_cluster[[1]]@match_cluster <- NA_character_

    expect_equal(matched_motif_2, matched_motif_no_cluster)
})
