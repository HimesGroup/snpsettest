##' Read in reference PLINK Bed file.
##'
##'
##' @param path path to the .bed file
##' @param check_monomorphic check monomorphic SNPs and remove them if exists
##' @param check_duplicates check duplicate SNPs by ID and genomic position, and
##'   remove them if exists
##' @param ... passed to [gaston::read.bed.matrix]
read_reference_bed <- function(path,
                               check_monomorphic = TRUE,
                               check_duplicates = TRUE,
                               ...) {

  ## to get basic stats for bed file
  op <- options(gaston.auto.set.stats = TRUE)

  path <- normalizePath(path)
  if (!file.exists(path)) {
    path <- paste0(path, ".bed")
    if (!file.exists(path)) {
      stop("Bed file not found.")
    }
  }
  ## remove extension
  path <- tools::file_path_sans_ext(path)
  x <- gaston::read.bed.matrix(basename = path, ...)

  if (check_monomorphic) {
    monomorphic <- gaston::test.snps(x, maf == 0)
    if (any(monomorphic)) {
      message("Removed ", sum(monomorphic), " monomoprhic SNPs.")
      x <- x[, !monomorphic]
    }
  }

  if (check_duplicates) {
    n1 <- ncol(x)
    x <- gaston::SNP.rm.duplicates(x, "id")
    n2 <- ncol(x)
    if ((n1 - n2) > 0) {
      message("Removed ", pretty_num((n2 - n1)), "duplicate SNPs by ID.")
    }
    x <- gaston::SNP.rm.duplicates(x, "chr:pos:alleles")
    n3 <- ncol(x)
    if ((n3 - n2) > 0) {
      message(
        "Removed ", pretty_num((n3 - n3)), "duplicate SNPs by genomic position."
      )
    }
  }

  ## Z-standardize genotype matrix
  ## if missing values exist in genotypes, we can impute them by;
  ## replacing NA with 0 -> imputing missing genotypes by the mean dosage
  gaston::standardize(x) <- "mu_sigma"
  message("Genotype matrix is standardized by ",
          "means and genotypic standard deviations.")

  ## reset options
  options(op)

  ## return bed.matrix
  x
}
