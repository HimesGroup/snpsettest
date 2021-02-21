##' Set-based association tests
##'
##' Perform set-based association tests between multiple sets of SNPs and a
##' phenotype using GWAS summary statistics.
##' @param hsumstats A data frame processed by \code{\link{harmonize_sumstats}}.
##' @param x A `bed.matrix` object created from the reference data.
##' @param snp_sets A named list where each index represents a separate set of
##'   SNPs.
##' @param method A method to compute a set-level p value.
##' @return A data.table with columns: "set.id", "p", "n.snp", "n.snp.clumped",
##'   "top.snp.id" and "top.snp.p"
##' - set.id = a name of SNP set
##' - p = a set-level p value
##' - n.snp = the number of SNPs in an intersection of set input and the
##' reference data
##' - n.snp.clumped = the number of SNPs in a set after LD clumping
##' - top.snp.id = SNP ID with the smallest p-value within a set of SNPs
##' - top.snp.p = The smallest p-value within a set of SNPs
##' @examples
##' ## Load GWAS summary data
##' data(exGWAS)
##'
##' ## Load reference genotype data
##' bfile <- system.file("extdata", "example.bed", package = "snpsettest")
##' x <- read_reference_bed(path = bfile)
##'
##' ## GWAS harmonization with reference data
##' hsumstats <- harmonize_sumstats(exGWAS, x)
##'
##' ## Load gene information data
##' data(gene.curated.GRCh37)
##'
##' ## Map SNPs to genes
##' snp_sets <- map_snp_to_gene(hsumstats, gene.curated.GRCh37)
##'
##' ## Perform set-based (gene-based) association tests
##' \dontrun{
##' snpset_test(hsumstats, x, snp_sets$sets)
##' }
##' @export
##' @importFrom stats qchisq
snpset_test <- function(hsumstats, x, snp_sets,
                        method = c("davies", "saddle")) {

  hsumstats_name <- deparse(substitute(hsumstats))

  is_df(hsumstats)
  is_bed_matrix(x)
  has_columns(hsumstats, c("id", "chr", "pos", "A1", "A2", "p"))
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
        gaston:: standardize(x) <- "mu_sigma"
      }
    )
  }

  ## Compute Chi-square stat from p
  hsumstats[, chisq := qchisq(p, df = 1, lower.tail = FALSE)]

  ## Perform set tests
  message("Starting set-based association tests...\n-----\n")
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
                     pd_tol = 1e-7, method = c("davies", "saddle")) {

  ## Assertion is not necessary since this function is only used internally.
  ## method <- match.arg(method)

  snp_ind <- match_cpp(snp_set, hsumstats$id)
  if (length(snp_ind) == 0L) {
    message("+++ Skipping ", set_id, ": no variants to test +++")
    return(NULL)
  }

  set_df <- hsumstats[snp_ind, ]
  set_p_min_ind <- which.min(set_df$p)
  top_snp_id <- set_df$id[set_p_min_ind]
  top_snp_p <- set_df$p[set_p_min_ind]

  message(
    "+++ Testing: ", set_id, " with ", pretty_num(nrow(set_df)),
    " SNPs +++"
  )

  ## Use set_df$id instead of snp_set for testing; there could be redundant SNPs
  ## in snp_set where they are not found in hsumstats but possibly in ref data.
  cor_ind <- match_cpp(set_df$id, x@snps$id)

  if (length(cor_ind) == 1L) {

    ## If only one SNP in a set, the set-based association p would be top_snp_p
    p <- top_snp_p

  } else {

    ## Extract genotype matrix; when missing_in_geno = `TRUE`, a returned matrix
    ## is imputed by mean dosage.
    ## Although correlation could be computed with pairwise deletion of missing
    ## values, it often causes negative eigenvalues so we prefer imputation.
    geno <- gaston::as.matrix(x[, cor_ind])
    if (missing_in_geno) {
      geno[is.na(geno)] <- 0 # remember geno is the Z-standardized matrix
    }

    cor_mat <- cor_cpp(geno)

    ## Get eigenvalues
    ev <- eigen(cor_mat, symmetric = TRUE, only.values = TRUE)$values

    ## Replacing negative or "almost zero" eigen values with a tolerance
    ## (adapted from sfsmisc::posdefify).
    ev[ev < pd_tol] <- pd_tol

    ## Calculate a test-statistic
    t_obs <- sum(set_df$chisq)
    if (method == "davies") {
      p <- tryCatch(
        davies(t_obs, lambda = ev, lim = 1e6, acc = 1e-8)$Qq,
        warning = function(w) {
          message(
            "- Davies method failed to produce a meaningful results. ",
            "Use Saddlepoint approximation instead."
          )
          pchisqsum(t_obs, df = rep(1, length(ev)), a = ev, lower.tail = FALSE)
        }
      )
    } else {
      p <- pchisqsum(t_obs, df = rep(1, length(ev)), a = ev, lower.tail = FALSE)
    }

  }
  ## Print set-based association p
  message("- P: ", prettyNum(p))

  ## Return output
  data.frame(
    set.id = set_id, p = p, n.snp = length(snp_ind),
    top.snp.id = top_snp_id, top.snp.p = top_snp_p
  )

}
