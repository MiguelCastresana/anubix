repo_root <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), ".."))
local_library <- file.path(repo_root, ".Rlib")

if (dir.exists(local_library)) {
  .libPaths(c(local_library, .libPaths()))
}

setwd(repo_root)
pkgload::load_all(repo_root, export_all = FALSE, helpers = FALSE)
testthat::test_dir(file.path(repo_root, "tests", "testthat"), reporter = "summary")
