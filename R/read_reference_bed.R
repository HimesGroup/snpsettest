##' Read in reference PLINK Bed file.
##'
##'
##' @param path A path to the .bed file
##' @param ... Further arguments used in [gaston::read.bed.matrix]
##' @export
read_reference_bed <- function(path, ...) {

  ## To get basic stats for bed file
  op <- options(gaston.auto.set.stats = TRUE)

  path <- normalizePath(path)
  if (!file.exists(path)) {
    path <- paste0(path, ".bed")
    if (!file.exists(path)) {
      stop("Bed file not found.")
    }
  }
  ## Remove .bed extension
  path <- tools::file_path_sans_ext(path)
  x <- gaston::read.bed.matrix(basename = path, ...)

  ## Find monomorphic SNPs
  monomorphic <- x@sigma == 0 # due to possible heterozyguous mono (maf = 0.05)
  if (any(monomorphic)) {
    message("Removed ", sum(monomorphic), " monomoprhic SNPs.")
    x <- x[, !monomorphic]
  }

  ## Find duplicate SNP IDs
  id_dup <- gaston::SNP.duplicated(x, "id")

  ## Find duplicate SNPs by base-pair and allele codes
  position_allele_dup <- gaston::SNP.duplicated(x, "chr:pos:alleles")

  ## Remove duplicates
  remove_ind <- unique(c(id_dup, position_allele_dup))
  if (length(remove_ind) > 0L) {
    keep_ind <- setdiff(seq_len(ncol(x)), remove_ind)
    message("Removed ", length(remove_ind), " duplicate SNPs.")
    x <- x[, keep_ind]
  }

  ## Z-standardize genotype matrix
  ## If missing values exist in genotypes, we can impute them by;
  ## replacing NA with 0 -> imputing missing genotypes by the mean dosage
  gaston::standardize(x) <- "mu_sigma"
  message("Genotype matrix is standardized by ",
          "means and genotypic standard deviations.")

  ## Reset options
  options(op)

  ## Return bed.matrix
  x
}
