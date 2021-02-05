##' Match alleles
##'
##' Match alleles between summary statistics and SNP information. This function
##' is a modified version of [bigsnpr::snp_match()], which removed unnecessary
##' requirements for set-based tests and improved speed with [data.table] syntax.
##' @noRd
##' @importFrom data.table melt setorder fcase
##' @importFrom vctrs vec_in
snp_match <- function(sumstats, info_snp, check_strand_flip = TRUE) {

  is_df(sumstats)
  is_df(info_snp)
  has_columns(sumstats, c("chr", "pos", "a1", "a2"))

  message2("%s variants to be matched.", format(nrow(sumstats), big.mark = ","))

  ## narrow down the scope of matching by physical position
  sumstats <- sumstats[vec_in(sumstats[, c("chr", "pos")],
                              info_snp[, c("chr", "pos")]),
                       c("chr", "pos", "a1", "a2", "p")]

  if (nrow(sumstats) == 0L) {
    stop2("No variant has been matched.")
  }

  sumstats[,  `:=` (a1_R = a2, a2_R = a1)]

  if (check_strand_flip) {
    sumstats <- remove_ambiguous_snps(sumstats)
    sumstats[, `:=` (a1_F = flip_strand(a1), a2_F = flip_strand(a2),
                     a1_RF = flip_strand(a1_R), a2_RF = flip_strand(a2_R))]
  }

  sumstats <-  melt(sumstats, measure = patterns("^a1", "^a2"),
                    value.name = c("a1", "a2"))

  sumstats <- info_snp[sumstats, on = c("chr", "pos", "a1", "a2"), nomatch = 0]

  message2(
    "%s variants have been matched.", format(nrow(sumstats), big.mark = ",")
  )

  setorder(sumstats[, variable := NULL], chr, pos)

}

flip_strand <- function(allele) {
  chartr("ACTG", "TGAC", allele)
}

remove_ambiguous_snps <- function(sumstats) {
    ambiguous_snps <- data.table(a1 = c("A", "T", "C", "G"),
                                 a2 = c("T", "A", "G", "C"))
    is_ambiguous <- vec_in(sumstats[, c("a1", "a2")], ambiguous_snps)
    message2(
      "%s ambiguous SNPs have been removed.",
      format(sum(is_ambiguous), big.mark = ",")
    )
    sumstats[!is_ambiguous, ]
}
