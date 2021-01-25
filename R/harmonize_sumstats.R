##' Harmonizing GWAS summary to reference data
##'
##' Finds an intersection of variants between GWAS summary and reference data.
##'
##' Pre-processing of GWAS summary data is required because the sets of variants
##' available to a particular GWAS might be poorly matched to the variants in
##' reference data. SNP matching is performed using chromosome code, base-pair
##' position, and allele codes, while taking into account possible strand flips
##' and reverse reference alleles. For matched entries, the SNP IDs in summary
##' data are replaced with the ones in the reference data.
##'
##' @param sumstats A data frame with columns: "snp.id", "chr", "pos", "a0",
##'   "a1" and "p".
##' - snp.id = SNP ID (e.g., rs numbers)
##' - chr =  chromosome (must be integer)
##' - pos =  base-pair position (must be integer)
##' - a0 = the first allele
##' - a1 = the second allele
##' - p = p value of SNP
##'
##' It could have only two columns "snp.id" and "p" if `join_by_id_only = TRUE`.
##' @param bigsnpobj A `bigSNP` object created using the reference data.
##' @param remove_ambiguous_snps `TRUE` or `FALSE`. If `TRUE`, the output
##'   removes ambiguous A/T and G/C SNPs for which the strand is not obvious. If
##'   the GWAS genotype data ifself is used as the reference data, it would be
##'   safe to set `FALSE`. Default is `TRUE`.
##' @param join_by_id_only `TRUE` or `FALSE`. If `TRUE`, SNP matching will be
##'   performed by SNP IDs only. Default is `FALSE`.
##' @return A data frame with columns: "snp.id", "chr", "pos", "a0", "a1" and
##'   "p".
##' @export
##' @importFrom data.table setDT setnames rbindlist :=
harmonize_sumstats <- function(sumstats, bigsnpobj,
                               join_by_id_only = FALSE,
                               remove_ambiguous_snps = TRUE
                               ) {

  sumstats_label <- deparse(substitute(sumstats))

  is_df(sumstats)
  check_class(bigsnpobj, "bigSNP")
  is_tf(remove_ambiguous_snps, "remove_ambiguous_snps")
  is_tf(join_by_id_only, "join_by_id_only")

  setDT(sumstats, key = "chr")
  info_snp <- get_snp_info(bigsnpobj)

  if (join_by_id_only) {
    has_columns(sumstats, c("snp.id", "p"))
    if (anyDuplicated(sumstats$snp.id) > 0) {
      stop2(paste0("Duplicate SNP IDs are found in %s. ",
                   "If join_by_id_only = TRUE, IDs should be unique."),
            sumstats_label)
    }
    ## sumstats <- merge(info_snp, sumstats[, c("snp.id", "p")],
    ##                   by = "snp.id", all = FALSE)
    ## setorder(sumstats, chr)
    sumstats <- info_snp[sumstats[, c("snp.id", "p")],
                         on = "snp.id", nomatch = 0]
    stopifnot(anyDuplicated(sumstats) == 0L) ## safety check
    message2("%s IDs are found in reference data.\n",
             format(nrow(sumstats), big.mark = ","))
  } else {
    has_columns(sumstats, c("chr", "pos", "a0", "a1", "p"))
    if (anyDuplicated(sumstats[, c("chr", "pos", "a0", "a1")])) {
      stop2(paste0("Some variants in %s ",
                   "share the same bp coordinate and allele codes."),
            sumstats_label)
    }
  }

  ## bigsnpr::snp_match requires a 'beta' column but it's not necessary for
  ## snpset_test
  sumstats$beta <- NA_integer_

  ## req_cols <- c("chr", "pos", "a0", "a1", "p", "beta")
  ## sumstats <- split(sumstats[, req_cols, with = FALSE], by = "chr")
  sumstats <- split(
    sumstats[, c("chr", "pos", "a0", "a1", "p", "beta")],
    by = "chr"
  )

  sumstats <- rbindlist(Map(
    function(x, i) {
      message2("+++ processing: chromosome %s +++", i)
      matched <- tryCatch(bigsnpr::snp_match(
        x,
        info_snp = info_snp,
        strand_flip = remove_ambiguous_snps,
        match.min.prop = 0
        ), error = function(e) {
          message2("no variants have been matched.")
          NULL
        })
    },
    sumstats,
    names(sumstats)
  ))

  ## sumstats$chr <- as.integer(sumstats$chr) ## sumstats[, chr := as.integer(chr)]

  stopifnot(anyDuplicated(sumstats$snp.id) == 0L) # safety check

  message2("\n%s variants have been matched in total.",
          format(nrow(sumstats), big.mark = ","))

  tryCatch(sumstats[, c("snp.id", "chr", "pos", "a0", "a1", "p")],
           error = function(e) invisible())

}

