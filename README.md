# ANUBIX

ANUBIX is an R package for network-based pathway enrichment analysis. It estimates the expected crosstalk between a query gene set and pathway gene sets by random sampling, then evaluates enrichment with a beta-binomial model.

<p align="center">
  <img src="Figure.png" alt="ANUBIX overview" width="420" />
</p>

For background and methodology, see the [ANUBIX paper](https://pubmed.ncbi.nlm.nih.gov/32788619/) and the PathBIX web implementation at <https://pathbix.sbc.su.se/>.

## What The Package Provides

| Function | Purpose |
| --- | --- |
| `anubix_links()` | Precomputes network links from each network gene to each pathway. Run this before enrichment. |
| `anubix()` | Runs the ANUBIX network-enrichment test with degree-constrained sampling. |
| `anubix_transitivity()` | Adds gene-set transitivity as an additional constraint. |
| `anubix_clustering()` | Clusters a query gene set with Infomap before applying ANUBIX. |

## Installation

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("MiguelCastresana/anubix")
library(ANUBIX)
```

## Basic Workflow

```r
library(ANUBIX)

links <- anubix_links(
  network = example_anubix$network,
  pathways = example_anubix$pathway_set,
  cutoff = 0.8,
  network_type = "weighted"
)

result <- anubix(
  network = example_anubix$network,
  links_matrix = links,
  genesets = example_anubix$gene_set,
  pathways = example_anubix$pathway_set,
  cores = 2,
  cutoff = 0.8,
  sampling = 2000,
  network_type = "weighted"
)
```

Important notes:

- Networks are treated as undirected.
- `anubix_links()` should be run before enrichment.
- Weighted networks are supported by `anubix_links()` when a weight column is provided.
- `anubix()` is the main enrichment entry point for most analyses.

## Development And Tests

Install test dependencies into the repo-local library:

```bash
Rscript tools/install-test-deps.R
```

Run unit tests:

```bash
Rscript tools/run-tests.R
```

Run the optional FunCoup plus KEGG integration test:

```bash
Rscript tools/run-integration-tests.R
```

The integration test downloads a small FunCoup network and KEGG pathway mappings for `Methanocaldococcus jannaschii`, then verifies that ANUBIX recovers an expected pathway signal.

## Repository Layout

```text
.
├── R/                 # Package source
├── man/               # Function documentation
├── tests/testthat/    # Unit and optional integration tests
├── tools/             # Dependency, test, and diagnostic helpers
├── DESCRIPTION
└── README.md
```

## Contact

Miguel Castresana Aguirre  
[miguel.castresana.aguirre@ki.se](mailto:miguel.castresana.aguirre@ki.se)
