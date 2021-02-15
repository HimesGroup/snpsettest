##' Match alleles
##'
##' Match alleles between summary statistics and SNP information. This function
##' is a modified version of [bigsnpr::snp_match()], which removed unnecessary
##' requirements for set-based tests and improved speed with [data.table] syntax.
##' @noRd
snp_match <- function(sumstats, info_snp, check_strand_flip = TRUE) {

  ## comment out assertions since this function is used only internally
  ## is_df(sumstats)
  ## is_df(info_snp)
  ## has_columns(sumstats, c("chr", "pos", "A1", "A2"))

  message(pretty_num(nrow(sumstats)), " variants to be matched.")

  ## narrow down the scope of SNP matching by physical position
  ind <- sumstats[, .(chr, pos)][info_snp[, .(chr, pos)], on = .(chr, pos),
                                 which = TRUE]
  sumstats <- sumstats[ind, .(chr, pos, A1, A2, p)]

  if (nrow(sumstats) == 0L) {
    stop("No variant has been matched.", call. = FALSE)
  }

  ## add columns to check reference allele swap
  sumstats[,  `:=` (A1_S = A2, A2_S = A1)]

  ## prepare an optional strand flip check
  if (check_strand_flip) {
    sumstats <- remove_ambiguous_snps(sumstats)
    sumstats[, `:=` (A1_F = flip_strand(A1),
                     A2_F = flip_strand(A2),
                     A1_SF = flip_strand(A1_S),
                     A2_SF = flip_strand(A2_S))]
  }

  sumstats <-  melt(sumstats, measure = patterns("^A1", "^A2"),
                    value.name = c("A1", "A2"))

  ## inner_join
  sumstats <- info_snp[sumstats, on = .(chr, pos, A1, A2), nomatch = NULL]

  message(pretty_num(nrow(sumstats)), " variants have been matched.")

  ## drop unnecessary variable from data.table::melt
  sumstats[, variable := NULL]

  ## return with sort
  setorder(sumstats, chr, pos)

}

flip_strand <- function(allele) {
  chartr("ACTG", "TGAC", allele)
}

remove_ambiguous_snps <- function(sumstats) {
  ambiguous_snps <- data.table(A1 = c("A", "T", "C", "G"),
                               A2 = c("T", "A", "G", "C"))

  ## get matching indices
  is_ambiguous <- sumstats[, .(A1, A2)][ambiguous_snps, on = .(A1, A2),
                                        which = TRUE]

  message(pretty_num(length(is_ambiguous)),
          " ambiguous SNPs have been removed.")

  sumstats[!is_ambiguous, ]
}
