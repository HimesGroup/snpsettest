##' Read a PLINK bed file for reference data
##'
##' Create a `bed.matrix` object from a .bed file. The function expects
##' .fam and .bim files under the same directory. See [gaston::read.bed.matrix]
##' for more details.
##' @param path A path to the .bed file
##' @param ... Further arguments used in [gaston::read.bed.matrix]
##' @return A [gaston::bed.matrix] object with a Z-standardized genotype matrix
##' @examples
##' ## Get a path to the example .bed file
##' bfile <- system.file("extdata", "example.bed", package = "snpsettest")
##'
##' ## Read a .bed file
##' x <- read_reference_bed(path = bfile)
##' @export
read_reference_bed <- function(path, ...) {

  ## To always get basic stats for bed file
  op <- options(gaston.auto.set.stats = TRUE)
  on.exit(options(op))

  path <- path.expand(path)
  if (!file.exists(path)) {
    path <- paste0(path, ".bed")
    if (!file.exists(path)) {
      stop("Bed file not found.")
    }
  }
  ## Remove .bed extension
  path <- tools::file_path_sans_ext(path)
  x <- gaston::read.bed.matrix(basename = path, ...)

  ## Z-standardize genotype matrix
  ## If missing values exist in genotypes, we can impute them by;
  ## replacing NA with 0 -> imputing missing genotypes by the mean dosage
  gaston::standardize(x) <- "mu_sigma"

  ## Return bed.matrix
  message("Created a bed.matrix with ",
          pretty_num(nrow(x)), " individuals and ",
          pretty_num(ncol(x)), " markers.")
  invisible(x)
}
