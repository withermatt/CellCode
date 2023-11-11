# CellCode

startup.R is sourced from the other Rmd files and loads the necessary R packages and user-defined functions.

cds_preprocess was used to do the quality contorl cell and gene filtering to generate the working cds file (cds_processed).

cds_by_base took cds_processed, split it by base condition, and then re-reduced the dimensions and clustered.

These 2 cds files are then used in the larger, analysis file.
