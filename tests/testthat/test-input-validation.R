test_that("anubix validates missing required inputs", {
  expect_error(
    anubix(NULL, data.frame(), data.frame(), data.frame()),
    "Network file is missing",
    fixed = TRUE
  )
  expect_error(
    anubix(data.frame(a = "A", b = "B"), NULL, data.frame(), data.frame()),
    "Precomputed file for the links per gene is missing",
    fixed = TRUE
  )
  expect_error(
    anubix(data.frame(a = "A", b = "B"), data.frame(), NULL, data.frame()),
    "Query gene sets are missing",
    fixed = TRUE
  )
  expect_error(
    anubix(data.frame(a = "A", b = "B"), data.frame(), data.frame(), NULL),
    "Pathway file is missing",
    fixed = TRUE
  )
})

test_that("anubix_transitivity validates missing required inputs", {
  expect_error(
    anubix_transitivity(NULL, data.frame(), data.frame(), data.frame()),
    "Network file is missing",
    fixed = TRUE
  )
  expect_error(
    anubix_transitivity(data.frame(a = "A", b = "B"), NULL, data.frame(), data.frame()),
    "Precomputed file for the links per gene is missing",
    fixed = TRUE
  )
  expect_error(
    anubix_transitivity(data.frame(a = "A", b = "B"), data.frame(), NULL, data.frame()),
    "Query gene sets are missing",
    fixed = TRUE
  )
  expect_error(
    anubix_transitivity(data.frame(a = "A", b = "B"), data.frame(), data.frame(), NULL),
    "Pathway file is missing",
    fixed = TRUE
  )
})

test_that("anubix_clustering validates missing required inputs", {
  expect_error(
    anubix_clustering(NULL, data.frame(), data.frame(), data.frame()),
    "Network file is missing",
    fixed = TRUE
  )
  expect_error(
    anubix_clustering(data.frame(a = "A", b = "B"), NULL, data.frame(), data.frame()),
    "Precomputed file for the links per gene is missing",
    fixed = TRUE
  )
  expect_error(
    anubix_clustering(data.frame(a = "A", b = "B"), data.frame(), NULL, data.frame()),
    "Query gene sets are missing",
    fixed = TRUE
  )
  expect_error(
    anubix_clustering(data.frame(a = "A", b = "B"), data.frame(), data.frame(), NULL),
    "Pathway file is missing",
    fixed = TRUE
  )
})
