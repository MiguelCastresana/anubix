# **ANUBIX**

**ANUBIX** is a genome-wide network analysis tool for pathway enrichment analysis. It is based on random sampling to build the expected crosstalk distribution between a query gene set and a pathway. The statistical significance is then assessed using a beta-binomial distribution.

For a detailed explanation of **ANUBIX** and its applications, please refer to the [ANUBIX paper](https://pubmed.ncbi.nlm.nih.gov/32788619/) .

### **Important Notes**
- **`anubix_links()`** function needs to be run before any other operations.
- **ANUBIX** is designed specifically for processing **undirected** networks.
- We recommend using the newer **`anubix_constrained`** function instead of **`anubix`** to obtain more comprehensive results.

## **Getting Started**

### **Installation**

To install **ANUBIX** from GitHub, use the following R code:

```r
# Install devtools or remotes if not already installed
install.packages("devtools")  # or install.packages("remotes")

# Install the ANUBIX package from GitHub
devtools::install_github("MiguelCastresana/anubix")  # or remotes::install_github("MiguelCastresana/anubix")

# Load the package
library(anubix)

