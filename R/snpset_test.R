##' Set-based association tests
##'
##' Perform set-based association tests between multiple sets of SNPs and a
##' phenotype using GWAS summary statistics. If the function encounters missing
##' genotypes in the reference data, they will be imputed with genotype means.
##' @param hsumstats A data frame processed by \code{\link{harmonize_sumstats}}.
##' @param x A `bed.matrix` object created from the reference data.
##' @param snp_sets A named list where each index represents a separate set of
##'   SNPs.
##' @param method A method to compute a set-level p value. "saddle" uses
##'   Kuonen's saddlepoint approximation (1999) and "davies" uses the algorithm
##'   of Davies (1980). Default is "saddle".
##' @return A data.table with columns: "set.id", "pvalue", "n.snp", "top.snp.id"
##'   and "top.snp.pvalue"
##' - set.id = a name of SNP set
##' - tstat = a test statistic
##' - pvalue = a set-level p value
##' - n.snp = the number of SNPs used in a test
##' - top.snp.id = SNP ID with the smallest p-value within a set of SNPs
##' - top.snp.pvalue = The smallest p-value within a set of SNP
##' @references
##' Kuonen, D. Saddlepoint Approximations for Distributions of Quadratic Forms
##' in Normal Variables. Biometrika 86, 929–935 (1999).
##'
##' Davies, R. B. Algorithm AS 155: The Distribution of a Linear Combination of
##' Chi-Square Random Variables. Journal of the Royal Statistical Society.
##' Series C (Applied Statistics) 29, 323–333 (1980).
##' @examples
##' ## GWAS summary statistics
##' head(exGWAS)
##'
##' ## Load reference genotype data
##' bfile <- system.file("extdata", "example.bed", package = "snpsettest")
##' x <- read_reference_bed(path = bfile)
##'
##' ## GWAS harmonization with reference data
##' hsumstats <- harmonize_sumstats(exGWAS, x)
##'
##' ## Perform a set-based test with an arbitrary SNP set
##' snpset_test(hsumstats, x, list(test = c("SNP_880", "SNP_1533", "SNP_4189")))
##'
##' ## Gene information data
##' head(gene.curated.GRCh37)
##'
##' ## Map SNPs to genes
##' snp_sets <- map_snp_to_gene(hsumstats, gene.curated.GRCh37)
##'
##' ## Perform gene-based association tests
##' \dontrun{
##' out <- snpset_test(hsumstats, x, snp_sets$sets)
##' }
##'
##' @export
##' @importFrom stats qchisq
snpset_test <- function(hsumstats, x, snp_sets,
                        method = c("saddle", "davies")) {

  hsumstats_name <- deparse(substitute(hsumstats))

  is_df(hsumstats)
  is_bed_matrix(x)
  has_columns(hsumstats, c("id", "chr", "pos", "A1", "A2", "pvalue"))
  is_named_list(snp_sets)
  method <- match.arg(method)

  if (!inherits(hsumstats, "data.table")) {
    setDT(hsumstats)
  }
  ## necessary ?
  if (!inherits(x@snps, "data.table")) {
    setDT(x@snps)
  }

  message("-----\n",
          pretty_num(nrow(hsumstats)),
          " variants are found in ",
          hsumstats_name,
          ".")
  message(pretty_num(length(snp_sets)),
          " set-based association tests will be performed.")

  ## Shrink hsumstats by removing SNPs not in test sets
  ## May not be good for memory usage since it makes a copy of hsumstats
  ## Could also get a subset of bed.matrix but think more about it...
  snps_in_sets <- unique(unlist(snp_sets))
  hsumstats <- hsumstats[id %in% snps_in_sets, ]

  ## Check any SNPs in hsumstats have missing genotypes in bed.matrix
  id_ind <- match_cpp(hsumstats$id, x@snps$id)
  missing_in_geno <- any(gaston::test.snps(x[, id_ind], callrate < 1L))

  ## If missing values exist, z-standardize genotypes so that:
  ## replacing NA with 0 == imputing missing genotypes by the mean dosage
  ## Z-standardization was done with read_reference_bed but check once again.
  if (missing_in_geno & !x@standardize_mu_sigma) {
    tryCatch(
      gaston::standardize(x) <- "mu_sigma", # propagate to subset by columns
      error = function(e) {
        x <- gaston::set.stats(x)
        gaston::standardize(x) <- "mu_sigma"
      }
    )
  }

  ## Compute Chi-square stat from p
  hsumstats[, chisq := qchisq(pvalue, df = 1, lower.tail = FALSE)]

  ## Perform set tests
  message("Starting set-based association tests...\n-----")
  rbindlist(
    Map(
      function(snp_set, set_id) {
        set_test(hsumstats, x, snp_set, set_id,
                 missing_in_geno = missing_in_geno, method = method)
      },
      snp_sets,
      names(snp_sets)
    )
  )
}

set_test <- function(hsumstats, x, snp_set, set_id, missing_in_geno,
                     pd_tol = 1e-7, method = c("saddle", "davies")) {

  ## Assertion is not necessary since this function is only used internally.
  ## method <- match.arg(method)

  snp_ind <- match_cpp(snp_set, hsumstats$id)
  if (length(snp_ind) == 0L) {
    message("+++ Skipping ", set_id, ": no variants to test +++")
    return(NULL)
  }

  set_df <- hsumstats[snp_ind, ]
  set_p_min_ind <- which.min(set_df$pvalue)
  top_snp_id <- set_df$id[set_p_min_ind]
  top_snp_p <- set_df$pvalue[set_p_min_ind]

  ## Use set_df$id instead of snp_set for testing; there could be redundant SNPs
  ## in snp_set where they are not found in hsumstats but possibly in ref data.
  cor_ind <- match_cpp(set_df$id, x@snps$id)

  if (length(cor_ind) == 1L) {

    ## If only one SNP in a set, the set-based association p would be top_snp_p
    p <- top_snp_p
    t_obs <- set_df$chisq

  } else {

    ## Extract genotype matrix; when missing_in_geno = `TRUE`, a returned matrix
    ## is imputed by mean dosage.
    ## Although correlation could be computed with pairwise deletion of missing
    ## values, it often causes negative eigenvalues so we prefer imputation.
    geno <- gaston::as.matrix(x[, cor_ind])
    if (missing_in_geno) {
      geno[is.na(geno)] <- 0 # remember geno is the Z-standardized matrix
    }

    ## Get eigenvalues using eigendecomposition of correlation matrix
    ## For a genotype matrix with n << p, use SVD instead
    if (ncol(geno) < 5 * nrow(geno)) {
      ev <- get_ev_from_evd(geno)
    } else {
      ev <- tryCatch(
        get_ev_from_svd(geno),
        error = function(e) get_ev_from_evd(geno)
      )
    }

    ## Replacing negative or "almost zero" eigen values with a tolerance
    ## (adapted from sfsmisc::posdefify).
    ev[ev < pd_tol] <- pd_tol

    ## Calculate a test-statistic
    t_obs <- sum(set_df$chisq)
    if (method == "saddle") {
      p <- pchisqsum(t_obs, df = rep(1, length(ev)), a = ev, lower.tail = FALSE)
    } else {
      p <- tryCatch(
        davies(t_obs, lambda = ev, lim = 1e6, acc = 1e-8)$Qq,
        warning = function(w) {
          pchisqsum(t_obs, df = rep(1, length(ev)), a = ev, lower.tail = FALSE)
        }
      )
    }

  }
  message(set_id, ": nSNP = ", pretty_num(length(cor_ind)), ", P = ",
          format(p, digit = 3))

  ## Return output
  data.table(
    set.id = set_id, tstat = t_obs, pvalue = p, n.snp = length(snp_ind),
    top.snp.id = top_snp_id, top.snp.pvalue = top_snp_p
  )

}
