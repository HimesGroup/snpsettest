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
##' @param sumstats A data frame with columns: "id", "chr", "pos", "A1",
##'   "A2" and "p".
##' - id = SNP ID (e.g., rs numbers)
##' - chr =  chromosome (must be integer)
##' - pos =  base-pair position (must be integer)
##' - A1, A2 = allele codes (allele order is not important)
##' - p = p value of SNP
##'
##' It could have only two columns "id" and "p" if `match_by_id = TRUE`.
##' @param x A `bed.matrix` object created using the reference data.
##' @param match_by_id `TRUE` or `FALSE`. If `TRUE`, SNP matching will be
##'   performed by SNP IDs instead of base-pair position and allele codes.
##'   Default is `FALSE`.
##' @param check_strand_flip `TRUE` or `FALSE`. If `TRUE`, the function 1)
##'   removes ambiguous A/T and G/C SNPs for which the strand is not obvious,
##'   and 2) attempts to find additional matching entries by flipping allele
##'   codes (i.e., A->T, T->A, C->G, G->A). If the GWAS genotype data ifself is
##'   used as the reference data, it would be safe to set `FALSE`. With
##'   `join_by_id_only = TRUE`, it simply determines whether to remove ambiguous
##'   SNPs. Default is `TRUE`.

##' @return A data frame with columns: "id", "chr", "pos", "A1", "A2" and
##'   "p".
##' @export
harmonize_sumstats <- function(sumstats, x,
                               match_by_id = FALSE,
                               check_strand_flip = TRUE
                               ) {

  ## Save the names of data input for error message
  sumstats_name <- deparse(substitute(sumstats))

  ## Assert function arguments
  is_df(sumstats)
  is_bed_matrix(x)
  is_tf(check_strand_flip, "check_strand_flip")
  is_tf(match_by_id, "match_by_id")

  ## Make data.table object
  if (!inherits(sumstats, "data.table")) {
    setDT(sumstats)
  }
  if (!inherits(x@snps, "data.table")) {
    setDT(x@snps)
  }

  if (match_by_id) {
    ## Check mandatory columns for `ID` matching
    has_columns(sumstats, c("id", "p"))

    ## Input SNP IDs must be unique
    if (anyDuplicated(sumstats$id) > 0L) {
      stop(
        "Duplicate SNP IDs are found in ",
        "'", sumstats_name, "'. ",
        "If '", sumstats_name, "' ",
        "is being matched by `id`, SNP IDs must be unique.",
        call. = FALSE
      )
    }

    message(pretty_num(nrow(sumstats)), " variants to be matched.")

    ## Inner join
    sumstats <- x@snps[sumstats[, .(id, p)], on = .(id), nomatch = NULL]

    ## Quick safety check; not strictly necessary
    stopifnot(anyDuplicated(sumstats) == 0L)

    ## Remove ambiguous SNPs and sort
    if (check_strand_flip) {
      message(pretty_num(nrow(sumstats)), " SNP IDs are shared.")
      sumstats <- remove_ambiguous_snps(sumstats)
      setorder(sumstats, chr, pos) # order by chr & pos
    }

    message(pretty_num(nrow(sumstats)),
            " variants have been matched.")

  } else {
    ## Check mandatory columns for `chr:pos:A1:A2` matching
    has_columns(sumstats, c("chr", "pos", "A1", "A2", "p"))

    ## Input `chr:pos:A1:A2` combination must be unique
    if (anyDuplicated(sumstats[, .(chr, pos, A1, A2)])) {
      stop(
        "Some variants in ",
        "'", sumstats_name, "' ",
        "share the same physical position and allele codes. ",
        "Remove duplicate SNPs in ",
        "'", sumstats_name, "'.",
        call. = FALSE
      )
    }

    ## Matching taking into account ref allele swaps (optionally strand flip)
    sumstats <- snp_match(sumstats, x@snps, check_strand_flip)
  }

  ## Return
  sumstats[, .(id, chr, pos, A1, A2, p)]

}

