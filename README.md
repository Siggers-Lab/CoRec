# hTF_array
Design, analysis, and visualization software for the human transcription factor (hTF) protein-binding microarray

## Installation

To make a local copy of the repository, navigate to the directory you want it in, and run `git clone https://github.com/Siggers-Lab/hTF_array.git`.

Then start an RStudio session (using R 4.0.5 if possible), navigate to the new repository directory, and double-click the "hTFArrayAnalysis.Rproj" file to open the project. Then run the following commands to install the necessary dependencies and load the package.

```
BiocManager::install("universalmotif")
devtools::load_all()
```

NOTE: If you are installing this package somewhere other than the BU SCC, you may need to install other dependencies as well using `install.packages("package_name")` or `BiocManager::install("package_name")` as appropriate.

## Example

To run the full analysis pipeline on the provided example data, run the following R command (also provided in `run_example.R`). The arguments are explained below.

```
corec_motifs <- run_full_analysis(
    output_directory = "./example_output",
    matrix_directory = "./example_data/data_matrices",
    matrix_base_name = "hTF_v1_SUDHL4_14jan21",
    pbm_conditions_file = "./example_data/v1_a11_run1_exp.txt",
    pbm_conditions = NA, # Need only one of pbm_conditions/pbm_conditions_file
    annotation_file = "./example_data/hTF_v01_PBM_ANNOT.txt",
    reference_motifs_file = "./example_data/JASPAR2018_hTF_only.meme",
    run_tag = "v1_a11_run1"
)
```

|        Argument       | Description |
|-----------------------|-------------|
|    output_directory   | The path to the directory in which to save the output files |
|    matrix_directory   | The path to the directory containing the individual o1, o2, or, and br fluorescence matrices to import. |
|    matrix_base_name   | The base filename of the fluorescence matrices to import. |
|  pbm_conditions_file  | A file containing the names of the PBM conditions, one per line. These will be used as column names in the fluorescence and z-score matrices. Only one of pbm_conditions_file or pbm_conditions should be provided. |
|     pbm_conditions    | A character vector of the names of the PBM conditions. These will be used as column names in the fluorescence and z-score matrices. Only one of pbm_conditions_file or pbm_conditions should be provided. |
|    annotation_file    | A file containing PBM probe annotations. |
| reference_motifs_file | A MEME format file containing the reference TF motifs to compare the hTF array motifs to. |
|        run_tag        | An optional (but recommended) tag specifying the particular experiment these fluorescence matrices are from. It will be incorporated into the column names along with the PBM conditions. This option is useful for differentiating replicates when merging two matrices downstream. |

After running the full analysis pipeline (which should take about 10-15 minutes), the output directory will contain two directories of .rds files of the generated corecmotif objects, a matrix of the raw fluorescence values, a matrix of the fluorescence z-scores, and a MEME format file containing the corecmotifs that passed the score threshold.
