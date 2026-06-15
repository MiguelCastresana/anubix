Sys.setenv(ANUBIX_RUN_INTEGRATION = "true")
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "run-tests.R"))
