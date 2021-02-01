<!-- badges: start -->
[![R-CMD-check](https://github.com/HimesGroup/snpsettest/workflows/R-CMD-check/badge.svg)](https://github.com/HimesGroup/snpsettest/actions)
<!-- badges: end -->

# snpsettest

The goal of the **snpsettest** package is to provide simple tools that perform a
set-based association test (e.g., gene-based association test) using GWAS
summary statistics. **This package is currently in alpha testing.**

To install this package, 

```R
devtools::install_github("HimesGroup/snpsettest")
```

A set-based test in the **snpsettest** package uses the statistical test
described in [VEGAS](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2896770/)
(**ve**rsatile **g**ene-based **a**ssociation **s**tudy), which combines the
effects of a set of SNPs within a set accounting for linkage disequilibrium
between markers. This test can be seen as a special case of
[SKAT](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3135811/#app3) (**s**equence
**k**ernel **a**ssociation **t**est) in which weights to the variants are all of
the same (i.e., unweighted SKAT).

<!-- Unlike the VEGAS software that used a simulation-based approach to calculate -->
<!-- p-values, this package used more computationally efficient methods as described -->
<!-- in here (to be added). -->

## TODO

- [ ] Write a user manual
- [ ] Test edge cases
- [ ] Remove redundancy in package import
- [ ] Add user friendly console output
- [ ] Re-expose some useful functions to end users
