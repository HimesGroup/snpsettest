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
##' @param sumstats A data frame with columns: "snp.id", "chr", "pos", "a1",
##'   "a2" and "p".
##' - snp.id = SNP ID (e.g., rs numbers)
##' - chr =  chromosome (must be integer)
##' - pos =  base-pair position (must be integer)
##' - a1, a2 = allele codes (allele order is not important)
##' - p = p value of SNP
##'
##' It could have only two columns "snp.id" and "p" if `join_by_id_only = TRUE`.
##' @param bigsnpobj A `bigSNP` object created using the reference data.
##' @param check_strand_flip `TRUE` or `FALSE`. If `TRUE`, the function 1)
##'   removes ambiguous A/T and G/C SNPs for which the strand is not obvious,
##'   and 2) attempts to find additional matching entries by flipping allele
##'   codes (i.e., A->T, T->A, C->G, G->A). If the GWAS genotype data ifself is
##'   used as the reference data, it would be safe to set `FALSE`. With
##'   `join_by_id_only = TRUE`, it simply determines whether to remove ambiguous
##'   SNPs. Default is `TRUE`.
##' @param join_by_id_only `TRUE` or `FALSE`. If `TRUE`, SNP matching will be
##'   performed by SNP IDs only. Default is `FALSE`.
##' @return A data frame with columns: "snp.id", "chr", "pos", "a1", "a2" and
##'   "p".
##' @export
##' @importFrom data.table setDT setnames rbindlist :=
harmonize_sumstats <- function(sumstats, bigsnpobj,
                               join_by_id_only = FALSE,
                               check_strand_flip = TRUE
                               ) {

  sumstats_label <- deparse(substitute(sumstats))

  is_df(sumstats)
  check_class(bigsnpobj, "bigSNP")
  is_tf(check_strand_flip, "check_strand_flip")
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

    sumstats <- info_snp[sumstats[, c("snp.id", "p")],
                         on = "snp.id", nomatch = 0]
    stopifnot(anyDuplicated(sumstats) == 0L) ## safety check
    message2("%s IDs are found in reference data.",
             format(nrow(sumstats), big.mark = ","))

    if (check_strand_flip) {
      sumstats <-  setorder(remove_ambiguous_snps(sumstats), chr, pos)
    }

    message2("%s variants have been matched.",
             format(nrow(sumstats), big.mark = ","))

  } else {
    has_columns(sumstats, c("chr", "pos", "a1", "a2", "p"))
    if (anyDuplicated(sumstats[, c("chr", "pos", "a1", "a2")])) {
      stop2(paste0("Some variants in %s ",
                   "share the same physical position and allele codes. ",
                   "Remove duplicate SNPs first."),
            sumstats_label)
    }

    ###############################
    ## per chromosome processing ##
    ###############################
    ## sumstats <- split(
    ##   sumstats[, c("chr", "pos", "a1", "a2", "p")],
    ##   by = "chr"
    ## )
    ## sumstats <- rbindlist(Map(
    ##   function(x, i) {
    ##     message2("+++ Processing: chromosome %s +++", i)
    ##     tryCatch(snp_match(
    ##       x,
    ##       info_snp = info_snp,
    ##       check_strand_flip = check_strand_flip
    ##     ), error = function(e) {
    ##       message2("No variant has been matched.")
    ##       NULL
    ##     })
    ##   },
    ##   sumstats,
    ##   names(sumstats)
    ## ))

    sumstats <- tryCatch(
      snp_match(sumstats, info_snp,
                check_strand_flip = check_strand_flip),
      error = function(e) {
        message2("No variant has been matched.")
        NULL
      })
  }

  stopifnot(anyDuplicated(sumstats$snp.id) == 0L) # safety check

  message2("\n%s variants have been matched in total.",
          format(nrow(sumstats), big.mark = ","))

  tryCatch(sumstats[, c("snp.id", "chr", "pos", "a1", "a2", "p")],
           error = function(e) invisible())

}

