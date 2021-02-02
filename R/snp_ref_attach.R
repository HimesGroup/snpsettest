##' Attach a "bigSNP" reference set from backing files into R
##'
##' This is a wrapper function for `bigsnpr::snp_attach()` with additional
##' validations for \code{\link{snpset_test}}.
##'
##' It ensures chromosome codes and base-pair positions are integers. If there
##' are duplicate SNP IDs in the reference set, it attempts to assign an unique
##' ID with a pattern: chr_pos_a1_a2. The updated IDs are saved back to
##' `rdsfile`.
##'
##' @param rdsfile A path of ".rds" which stores the `bigSNP` object.
##' @return A `bigSNP` object.
##' @export
##' @importFrom vctrs vec_duplicate_detect
##' @seealso [bigsnpr::snp_attach()] for more details.
snp_ref_attach <- function(rdsfile) {

  ## convert a file path to canonical form
  rdsfile <- normalizePath(rdsfile)

  ## attach a .rds file with bigsnpr::snp_attach()
  bigsnpobj <- bigsnpr::snp_attach(rdsfile)

  ## ensure chromosome codes are integer
  if (!(is_integer_vector(bigsnpobj$map$chromosome))) {
    stop2(paste0(
      "Chromosome codes must be integer values. ",
      "You may use PLINK --output-chr option to specify",
      "a numerical coding scheme in your reference data.")
      )
  }

  ## ensure base-pair positions are integer
  if (!(is_integer_vector(bigsnpobj$map$physical.pos))) {
    stop2("Base-pair positions must be integer values.")
  }

  ## assign unique SNP IDs if there are duplicates
  if (anyDuplicated(bigsnpobj$map$marker.ID) > 0L) {
    warning2(paste0(
      "Some SNP IDs in reference data are not unique. ",
      "They will be replaced with new IDs 'CHR_POS_A1_A2'."
    ))
    dup_idx <- get_duplicate_indice(bigsnpobj$map$marker.ID)
    bigsnpobj$map$marker.ID[dup_idx] <- paste(
      bigsnpobj$map$chromosome[dup_idx],
      ## as.integer for preventing paste with e notation (e.g., 2e+07)
      as.integer(bigsnpobj$map$physical.pos[dup_idx]),
      bigsnpobj$map$allele2[dup_idx],
      bigsnpobj$map$allele1[dup_idx],
        sep = "_"
      )

    if (anyDuplicated(bigsnpobj$map$marker.ID)) {
      stop2(paste0(
        "Cannot assign unique IDs by base-pair coordinate and allele codes. ",
        "Please check the source file again to remove duplicate variants."))
    }
    bigsnpr::snp_save(bigsnpobj) # save the updated IDs
    message2("The updated IDs are saved back to %s.", rdsfile)
    bigsnpobj <- bigsnpr::snp_attach(rdsfile) # reload .rds file
  }
  bigsnpobj
}
