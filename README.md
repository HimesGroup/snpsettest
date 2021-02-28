
<!-- README.md is generated from README.Rmd. Please edit that file -->

# snpsettest

<!-- badges: start -->

[![R-CMD-check](https://github.com/HimesGroup/snpsettest/workflows/R-CMD-check/badge.svg)](https://github.com/HimesGroup/snpsettest/actions)
<!-- badges: end -->

The goal of **snpsettest** is to provide simple tools that perform a
set-based association test (e.g., gene-based association test) using
GWAS summary statistics. **This package is currently under
development.**

A set-based association test in the **snpsettest** package is based on
the statistical test described in
[VEGAS](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2896770/)
(**ve**rsatile **g**ene-based **a**ssociation **s**tudy), which combines
the effects of a set of SNPs accounting for linkage disequilibrium
between markers. This package uses a more efficient method than the
original VEGAS implementation for computing a set-level p-value,
significantly reducing computation time.

## Installation

To install this package,

``` r
devtools::install_github("HimesGroup/snpsettest")
```

## Getting started

This is a basic example which shows you how to perform gene-based
association tests using GWAS summary statistics in which sets of SNPs
are defined by genes.

### GWAS summary file

The **snpsettest** package requires SNP-level p-values to perform
gene-based association tests.

``` r
library(snpsettest)

# Load an example of GWAS summary file
data(exGWAS)
head(exGWAS, 3)
#>      id chr   pos A1 A2    pvalue
#> 1 SNP_0   1 50215  G  C 0.1969353
#> 2 SNP_2   1 50768  A  G 0.6620465
#> 3 SNP_3   1 50833  T  G 0.5822596
```

### Reference data

To infer the relationships among SNPs, the **snpsettest** package
requires a reference data set. The GWAS genotype data itself can be used
as the reference data (If the GWAS cohort is large, it is impractical to
use genotype data of all individuals. It would be sufficient to randomly
select 1,000 unrelated individuals for inferring pairwise LD
correlations among common SNPs). Otherwise, you could use publicly
available data, such as the 1000 Genomes (please see the companion
vignette for [processing the 1000 Genomes
data](reference_1000Genomes.html)). This package accepts PLINK 1 binary
files (.bed, .bim, .fam) as an input. We can use `read_reference_bed` to
read them into R.

``` r
# Path to .bed file
bfile <- system.file("extdata", "example.bed", package = "snpsettest")

# Read a .bed file using bed.matrix-class in gaston package
# Genotypes are retrieved on demand to manage large-scale genotype data
x <- read_reference_bed(bfile, verbose = FALSE)
```

### Harmonize GWAS summary to the reference data

Pre-processing of GWAS summary data is required because the sets of
variants available in a particular GWAS might be poorly matched to the
variants in reference data. SNP matching can be performed using
`harmonize_sumstats` either 1) by SNP ID or 2) by chromosome code,
base-pair position, and allele codes, while taking into account
reference allele swap and possible strand flips.

``` r
# Harmonize by SNP IDs
hsumstats1 <- harmonize_sumstats(exGWAS, x)
#> -----
#> Checking the reference data for harmonization...
#> Found 0 monomoprhic SNPs in the reference data.
#> Found 0 duplicate SNP IDs in the reference data.
#> Excluded 0 SNPs from the harmonization.
#> -----
#> Checking the GWAS summary statistics...
#> 2,753 variants to be matched.
#> 2,630 variants have been matched.

# Harmonize by genomic position and allele codes
# Reference allele swap will be taken into account (e.g., A/C match C/A)
hsumstats2 <- harmonize_sumstats(exGWAS, x, match_by_id = FALSE)
#> -----
#> Checking the reference data for harmonization...
#> Found 0 monomoprhic SNPs in the reference data.
#> Found 0 duplicate SNPs in the reference data by genomic position and alleles codes.
#> Excluded 0 SNPs from the harmonization.
#> -----
#> Checking the GWAS summary statistics...
#> 2,753 variants to be matched.
#> 2,618 variants have been matched.

# Check matching entries by flipping allele codes (e.g., A/C match T/G)
# Ambiguous SNPs will be excluded from harmonization
hsumstats3 <- harmonize_sumstats(exGWAS, x, match_by_id = FALSE, check_strand_flip = TRUE)
#> -----
#> Checking the reference data for harmonization...
#> Found 0 monomoprhic SNPs in the reference data.
#> Found 0 duplicate SNPs in the reference data by genomic position and alleles codes.
#> Excluded 0 SNPs from the harmonization.
#> -----
#> Checking the GWAS summary statistics...
#> 2,753 variants to be matched.
#> 835 ambiguous SNPs have been removed.
#> 1,795 variants have been matched.
```

### Map SNPs to genes

To perform gene-based association tests, it is necessary to annotate
SNPs onto their neighboring genes. Mapping SNPs to genes (or genomic
regions) can be achieved by `map_snp_to_genes` with gene start/end
information.

``` r
# Load gene information
data(gene.curated.GRCh37) # extracted from GENCODE release 19
head(gene.curated.GRCh37, 3)
#>             gene.id chr  start    end strand  gene.name      gene.type
#> 1 ENSG00000186092.4   1  69091  70008      +      OR4F5 protein_coding
#> 2 ENSG00000237683.5   1 134901 139379      - AL627309.1 protein_coding
#> 3 ENSG00000235249.1   1 367640 368634      +     OR4F29 protein_coding

# Map SNPs to genes
snp_sets <- map_snp_to_gene(hsumstats1, gene.curated.GRCh37)
str(snp_sets$sets[1:5])
#> List of 5
#>  $ ENSG00000186092.4: chr [1:110] "SNP_0" "SNP_2" "SNP_3" "SNP_4" ...
#>  $ ENSG00000237683.5: chr [1:109] "SNP_317" "SNP_320" "SNP_321" "SNP_323" ...
#>  $ ENSG00000235249.1: chr [1:95] "SNP_1283" "SNP_1285" "SNP_1287" "SNP_1288" ...
#>  $ ENSG00000185097.2: chr [1:96] "SNP_2392" "SNP_2396" "SNP_2397" "SNP_2398" ...
#>  $ ENSG00000187634.6: chr [1:135] "SNP_3455" "SNP_3456" "SNP_3458" "SNP_3459" ...

# Allows a certain kb window before/after the gene to be included for SNP mapping
snp_sets_50kb <- map_snp_to_gene(
  hsumstats1, gene.curated.GRCh37, 
  extend_start = 50, extend_end = 50 # default is 20kb
)
```

### Perform gene-based association tests

Once we have SNP sets for genes, `snpset_test` can be used to perform
gene-based association tests.

``` r
# Perform gene-based association tests for the first 5 genes
res <- snpset_test(hsumstats1, x, snp_sets$sets[1:5])
#> -----
#> 2,630 variants are found in hsumstats1.
#> 5 set-based association tests will be performed.
#> Starting set-based association tests...
#> -----
#> ENSG00000186092.4: nSNP = 110, P = 0.042
#> ENSG00000237683.5: nSNP = 109, P = 0.00936
#> ENSG00000235249.1: nSNP = 95, P = 0.182
#> ENSG00000185097.2: nSNP = 96, P = 0.122
#> ENSG00000187634.6: nSNP = 135, P = 0.0103

# Show output
res
#>               set.id    tstat      pvalue n.snp top.snp.id top.snp.pvalue
#> 1: ENSG00000186092.4 141.7800 0.042027775   110     SNP_78   0.0009143436
#> 2: ENSG00000237683.5 154.2858 0.009362739   109    SNP_363   0.0006419257
#> 3: ENSG00000235249.1 109.0270 0.182400780    95   SNP_1311   0.0047610286
#> 4: ENSG00000185097.2 114.7301 0.122042173    96   SNP_2458   0.0034444534
#> 5: ENSG00000187634.6 185.7576 0.010306441   135   SNP_3601   0.0003350840
```
