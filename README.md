
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
association tests using GWAS summary statistics:

### GWAS summary file

``` r
library(snpsettest)

# Load an example of GWAS summary file
# snpsettest requires id, chr, pos, A1, A2, and p columns for GWAS summary file
data(exGWAS)
head(exGWAS)
#>      id chr   pos A1 A2          p
#> 1 SNP_0   1 50149  A  C 0.73013485
#> 2 SNP_2   1 50818  T  C 0.87234011
#> 3 SNP_3   1 51094  C  T 0.52959587
#> 4 SNP_5   1 51476  A  C 0.86017713
#> 5 SNP_6   1 51820  C  A 0.81271456
#> 6 SNP_7   1 51880  C  G 0.06055741
```

### Reference data

To infer the relationships among SNPs, the **snpsettest** package
requires a reference data set. The GWAS genotype data itself can be used
as the reference data. Otherwise, you could use publicly available data,
such as the 1000 Genome. This package accepts PLINK 1 binary files
(.bed, .bim, .fam) as an input.

``` r
# Path to .bed file
bfile <- system.file("extdata", "example.bed", package = "snpsettest")

# Read a .bed file using bed.matrix-class in gaston package
# Genotypes are retrieved on demand to manage large-scale genotype data
x <- read_reference_bed(bfile)
#> Genotype matrix was standardized by means and genotypic standard deviations.
```

### Harmonize GWAS summary to the reference data

Pre-processing of GWAS summary data is required because the sets of
variants available in a particular GWAS might be poorly matched to the
variants in reference data. SNP matching can be performed either 1) by
SNP ID or 2) by chromosome code, base-pair position, and allele codes,
while taking into account possible strand flips and reference allele
swap.

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
#> Found 0 duplicate SNPs in the reference data by base-pair position and alleles codes.
#> Excluded 0 SNPs from the harmonization.
#> -----
#> Checking the GWAS summary statistics...
#> 2,753 variants to be matched.
#> 2,625 variants have been matched.

# Check matching entries by flipping allele codes (e.g., A/C match T/G)
# Ambiguous SNPs will be excluded from harmonization
hsumstats3 <- harmonize_sumstats(exGWAS, x, match_by_id = FALSE, check_strand_flip = TRUE)
#> -----
#> Checking the reference data for harmonization...
#> Found 0 monomoprhic SNPs in the reference data.
#> Found 0 duplicate SNPs in the reference data by base-pair position and alleles codes.
#> Excluded 0 SNPs from the harmonization.
#> -----
#> Checking the GWAS summary statistics...
#> 2,753 variants to be matched.
#> 859 ambiguous SNPs have been removed.
#> 1,771 variants have been matched.
```

### Map SNPs to genes

To perform gene-based association tests, it is necessary to annotate
SNPs onto their neighboring genes.

``` r
# Load gene information
data(gene.curated.GRCh37) # extracted from gencode release 19

# Map SNPs to genes
snp_sets <- map_snp_to_gene(hsumstats1, gene.curated.GRCh37)
str(snp_sets$sets[1:5])
#> List of 5
#>  $ ENSG00000186092.4: chr [1:92] "SNP_0" "SNP_2" "SNP_3" "SNP_5" ...
#>  $ ENSG00000237683.5: chr [1:112] "SNP_307" "SNP_308" "SNP_310" "SNP_311" ...
#>  $ ENSG00000235249.1: chr [1:99] "SNP_1303" "SNP_1305" "SNP_1307" "SNP_1309" ...
#>  $ ENSG00000185097.2: chr [1:89] "SNP_2392" "SNP_2394" "SNP_2395" "SNP_2397" ...
#>  $ ENSG00000187634.6: chr [1:116] "SNP_3440" "SNP_3442" "SNP_3448" "SNP_3449" ...

# Allows a certain kb window before/after the gene to be included for SNP mapping
snp_sets_50kb <- map_snp_to_gene(
  hsumstats1, gene.curated.GRCh37, 
  extend_start = 50, extend_end = 50 # default is 20kb
)
```

### Perform gene-based association tests

``` r
# Perform gene-based association tests for the first 5 genes
res <- snpset_test(hsumstats1, x, snp_sets$sets[1:5])
#> 
#> -----
#> 2,630 variants are found in hsumstats1.
#> 5 set-based association tests will be performed.
#> Starting set-based association tests...
#> -----
#> +++ Testing: ENSG00000186092.4 with 92 SNPs +++
#> - P: 0.006626555
#> +++ Testing: ENSG00000237683.5 with 112 SNPs +++
#> - P: 0.002337155
#> +++ Testing: ENSG00000235249.1 with 99 SNPs +++
#> - P: 0.329092
#> +++ Testing: ENSG00000185097.2 with 89 SNPs +++
#> - P: 0.1582166
#> +++ Testing: ENSG00000187634.6 with 116 SNPs +++
#> - P: 0.01333011

# Show output
res
#>               set.id           p n.snp top.snp.id    top.snp.p
#> 1: ENSG00000186092.4 0.006626555    92    SNP_162 0.0009888022
#> 2: ENSG00000237683.5 0.002337155   112    SNP_436 0.0029226968
#> 3: ENSG00000235249.1 0.329092035    99   SNP_1440 0.0003399570
#> 4: ENSG00000185097.2 0.158216569    89   SNP_2397 0.0170853467
#> 5: ENSG00000187634.6 0.013330114   116   SNP_3524 0.0032527617
```
