##' Set-based association tests.
##'
##' Perform set-based association tests between multiple sets of SNPs and a
##' phenotype using GWAS summary statistics.
##' @param hsumstats A data frame processed by \code{\link{harmonize_sumstats}}.
##' @param bigsnpobj A `bigSNP` object created from the reference data.
##' @param snp_sets A named list where each index represents a separate set of
##'   SNPs.
##' @param missing_in_geno `TRUE` or `FALSE`. It determines if genotypes in the
##'   reference data have missing values. If `TRUE`, missing genotypes are
##'   imputed by the mode of each SNP. Default is 'TRUE'.
##' @param thr_rs A LD clumping threshold over the squared correlation between
##'   two SNPs. Default is 0.8.
##' @return A data.table with columns: "set.id", "p", "n.snp", "n.snp.clumped",
##'   "top.snp.id" and "top.snp.p"
##' - set.id = a name of SNP set
##' - p = a set-based p value
##' - n.snp = the number of SNPs in a set of SNPs
##' - n.snp.clumped = the number of SNPs in a set of SNPs after LD clumping
##' - top.snp.id = SNP ID with the smallest p-value within a set of SNPs
##' - top.snp.p = The smallest p-value within a set of SNPs
##' @export
##' @importFrom data.table setDT rbindlist :=
##' @importFrom Matrix nearPD
##' @importFrom stats qchisq
snpset_test <- function(hsumstats, bigsnpobj, snp_sets,
                        missing_in_geno = TRUE,
                        thr_rs = 0.8) {

  is_df(hsumstats)
  has_columns(hsumstats, c("snp.id", "chr", "pos", "a0", "a1", "p"))
  check_class(bigsnpobj, "bigSNP")
  is_named_list(snp_sets)
  is_number_between(thr_rs, 0L, 1L, "thr_rs")

  if (inherits(hsumstats, "data.table")) {
    setDT(hsumstats)
  }

  hsumstats[, chisq := qchisq(p, df = 1, lower.tail = FALSE)]

  info_snp <- get_snp_info(bigsnpobj)
  info_snp[hsumstats, on = intersect(names(info_snp), names(hsumstats)),
           `:=`(p = i.p, chisq = i.chisq)]

  G <- bigsnpobj$genotypes
  if (missing_in_geno) {
    G_noNA <- bigsnpr::snp_fastImputeSimple(G, method = "mean2")
  } else {
    G_noNA <- G # same mem address
  }

  ## temporary safety check
  stopifnot(identical(ncol(G), nrow(info_snp)))
  stopifnot(identical(info_snp$snp.id, bigsnpobj$map$marker.ID))
  stopifnot(identical(info_snp$chr, bigsnpobj$map$chromosome))
  stopifnot(identical(info_snp$pos, bigsnpobj$map$physical.pos))

  ## perform set tests
  rbindlist(
    Map(
      function(x, i) {
      set_test(info_snp, G, G_noNA, snp_set = x, set_name = i, thr_rs = thr_rs)
    },
    snp_sets,
    names(snp_sets)
    )
  )
}

set_test <- function(info_snp, G, G_noNA, snp_set, set_name, thr_rs) {

  message2("testing: %s", set_name)
  snp_ind <- which(info_snp$snp.id %in% snp_set)
  top_snp_p <- min(info_snp$p[snp_ind])
  top_snp_id <- snp_set[which.min(info_snp$p[snp_ind])]

  if (thr_rs < 1L) {
    cor_ind <- clumping(info_snp, G_noNA, snp_set, thr_rs)
  } else {
    cor_ind <- which(info_snp$snp.id %in% snp_set)
  }

  t_obs <- sum(info_snp$chisq[cor_ind])
  cor_mat <- bigsnpr::snp_cor(G, ind.col = cor_ind, size = Inf)
  ev <- eigen(cor_mat, symmetric = TRUE, only.values = TRUE)$values

  if(!is_positive_definite(ev)) {
    ev <- nearPD(cor_mat, corr = TRUE, only.values = TRUE)
  }

  out <- davies(t_obs, lambda = ev, h = rep(1, length(cor_ind)),
                lim = 1e7, acc = 1e-9)
  data.table(set.id = set_name, p = out$Qq, n.snp = length(snp_ind),
             n.snp.clumped = length(cor_ind),
             top.snp.id = top_snp_id, top.snp.p = top_snp_p)
}

clumping <- function(info_snp, G_noNA, snp_set,
                     thr_r2 = 0.8) {
  check_class(G_noNA, "FBM.code256")

  ind_exclude <- clumping_exclude_indices(
    info_snp$snp.id, snp_set
  )
  bigsnpr::snp_clumping(
    G = G_noNA,
    thr.r2 = thr_r2,
    S = info_snp$chisq,
    infos.chr = info_snp$chr,
    infos.pos = info_snp$pos,
    exclude = ind_exclude,
  )
}

clumping_exclude_indices <- function(superset, set) {
  which(!(superset %in% set))
}
