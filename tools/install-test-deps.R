packages <- c(
  "BiocManager",
  "dplyr",
  "nloptr",
  "optimx",
  "TailRank",
  "igraph",
  "purrr",
  "tibble",
  "testthat"
)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
repo_root <- normalizePath(file.path(dirname(script_path), ".."))
library_path <- file.path(repo_root, ".Rlib")
dir.create(library_path, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(library_path, .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = library_path, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("Biobase", quietly = TRUE)) {
  BiocManager::install("Biobase", lib = library_path, ask = FALSE, update = FALSE)
}

archive_packages <- c(
  "https://cran.r-project.org/src/contrib/Archive/oompaBase/oompaBase_3.2.9.tar.gz",
  "https://cran.r-project.org/src/contrib/Archive/oompaData/oompaData_3.1.4.tar.gz"
)

install.packages(archive_packages, lib = library_path, repos = NULL, type = "source")
install.packages(packages, lib = library_path, repos = "https://cloud.r-project.org")
