##' Get SNP information from a `bigSNP` object.
##'
##' @param bigsnpobj A `bigSNP` object created from the reference data.
##' @importFrom data.table data.table setDT setnames
##' @noRd
get_snp_info <- function(bigsnpobj) {
  ## for global variable problem due to non-standard evaluation
  ## pre-assign NULL at the top of function.
  ## snp.id <- chr <- pos <- a0 <- a1 <- NULL
  check_class(bigsnpobj, "bigSNP")
  info_snp <- bigsnpobj$map
  setDT(info_snp)
  setnames(info_snp,
    old = c("chromosome", "marker.ID", "physical.pos", "allele1", "allele2"),
    new = c("chr", "snp.id", "pos", "a1", "a0")
  )
  info_snp[, list(snp.id, chr, pos, a0, a1)]
}
