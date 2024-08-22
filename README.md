## ANUBIX ##
ANUBIX is a genome-wide network analysis tool for pathway enrichment analysis. It is based on random sampling to build the expected crosstalk distribution between a query gene set and a pathway. The statistical significance is then assessed using a beta-binomial distribution.

## Important Notes 
anubix_links() function needs to be run before any other operations.
ANUBIX is designed specifically for processing undirected networks.
We recommend using the newer anubix_constrained function instead of anubix to obtain improved results.

## Instalation ##

require(devtools)
devtools::install_github("MiguelCastresana/anubix")

Once it is installed follow the manual of the package ANUBIX (ANUBIX_MANUAL.pdf)
