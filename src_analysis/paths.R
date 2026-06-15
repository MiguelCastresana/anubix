# Shared paths for ANUBIX analysis scripts.

script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[1]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

benchmark_root <- Sys.getenv(
  "ANUBIX_BENCHMARK_DIR",
  unset = file.path(repo_root, "anubix_benchmark")
)

output_root <- Sys.getenv(
  "ANUBIX_OUTPUT_DIR",
  unset = file.path(repo_root, "results")
)

benchmark_file <- function(...) {
  file.path(benchmark_root, ...)
}

output_file <- function(...) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  file.path(output_root, ...)
}
