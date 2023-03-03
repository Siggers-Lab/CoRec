# CoRec
CoRec (short for Cofactor Recruitment) is a protein binding microarray (PBM) based approach for profiling transcription cofactor (COF) recruitment to particular DNA sequences. This package provides tools for analyzing and visualizing the data from CoRec experiments.

## Installation

The `CoRec` package requires a minimum R version of 4.1. It also requires an installation of the [MEME suite](https://meme-suite.org/meme/) (see ["Detecting the MEME Suite"](#detecting-the-meme-suite) for more information).

```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("Siggers-Lab/CoRec", build_vignettes = TRUE)
```

## Getting started

The `Introduction to CoRec` vignette walks through an example CoRec analysis.

```
library(CoRec)
vignette("CoRec")
```

## Detecting the MEME Suite

The `CoRec` package depends on the `memes` package for running motif comparisons. `memes` is an R wrapper for the [MEME suite](https://meme-suite.org/meme/), and it requires that the MEME suite be installed and findable. See the `memes` package ["Install MEME" vignette](https://snystrom.github.io/memes-manual/articles/install_guide.html) for more information.

For BU SCC users, the most convenient way to ensure that `memes` can find the MEME suite installation is to create or modify a `.Renviron` file in your home directory and add the following line:

```
MEME_BIN=/share/pkg.7/meme/5.3.3/install/bin/
```

Another option is to provide the path every time you run `process_corecmotifs()` or `find_match()`, e.g.:

```
process_corecmotifs(
    meme_path = "/share/pkg.7/meme/5.3.3/install/bin/",
    ...
)
```
