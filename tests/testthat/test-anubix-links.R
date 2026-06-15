test_that("anubix_links counts links in unweighted networks", {
  network <- data.frame(
    gene_a = c("A", "A", "B", "C", "D"),
    gene_b = c("B", "C", "D", "D", "E")
  )
  pathways <- data.frame(
    gene = c("B", "D", "E", "Z"),
    pathway = c("P1", "P1", "P2", "P1")
  )

  links <- anubix_links(
    network = network,
    pathways = pathways,
    network_type = "unweighted"
  )

  expect_s3_class(links, "data.frame")
  expect_named(links, c("P1", "P2"))
  expect_equal(rownames(links), c("A", "B", "C", "D", "E"))
  expect_equal(links["A", "P1"], 1)
  expect_equal(links["B", "P1"], 1)
  expect_equal(links["D", "P2"], 1)
  expect_equal(links["E", "P1"], 1)
})

test_that("anubix_links sums edge weights after applying the cutoff", {
  network <- data.frame(
    gene_a = c("A", "A", "B", "C", "D"),
    gene_b = c("B", "C", "D", "D", "E"),
    weight = c(0.90, 0.70, 0.85, 0.95, 0.90)
  )
  pathways <- data.frame(
    gene = c("B", "D", "E"),
    pathway = c("P1", "P1", "P2")
  )

  links <- anubix_links(
    network = network,
    pathways = pathways,
    cutoff = 0.8,
    network_type = "weighted"
  )

  expect_equal(links["A", "P1"], 0.90)
  expect_equal(links["C", "P1"], 0.95)
  expect_equal(links["D", "P1"], 0.85)
  expect_equal(links["D", "P2"], 0.90)
})

test_that("anubix_links validates required inputs", {
  pathways <- data.frame(gene = "A", pathway = "P1")
  network <- data.frame(gene_a = "A", gene_b = "B")

  expect_error(
    anubix_links(NULL, pathways),
    "Network file is missing",
    fixed = TRUE
  )
  expect_error(
    anubix_links(network, NULL),
    "Pathways file is missing",
    fixed = TRUE
  )
  expect_error(
    anubix_links(network, pathways, cutoff = "0.8"),
    "Link confidence cutoff is not in a proper format",
    fixed = TRUE
  )
  expect_error(
    anubix_links(network, pathways, network_type = "weighted"),
    "A weighted network must contain a third column",
    fixed = TRUE
  )
  expect_error(
    anubix_links(network, pathways, network_type = "directed"),
    "network_type must be either",
    fixed = TRUE
  )
})
