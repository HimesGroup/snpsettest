##' Map SNPs to genomic regions (e.g., genes)
##'
##' Map SNPs to genomic regions to extract sets of SNPs for set-based tests. A
##' typical usage would be to annotate SNPs onto their neighboring genes to
##' perform gene-based tests.
##'
##' @param info_snp A data frame with columns: "id", "chr", and "pos".
##' - id = an unique identifier (e.g., rs numbers)
##' - chr = chromosome
##' - pos = base-pair position
##' @param info_region A data frame with columns: "region.id", "chr", "start",
##'   and "end".
##' - region.id = an unique identifier (e.g., gene names)
##' - chr = chromosome (should use the same chromosome scheme in `info_snp`)
##' - start = start position of genomic region
##' - end =  end position of genomic region
##' @param extend_region_start A single non-negative integer, allowing for a
##'   certain kb window before the region to be included. Default is 20 (= 20kb).
##' @param extend_region_end A single non-negative integer, allowing for a
##'   certain kb window after the region to be included. Default is 20 (= 20kb).
##' @param only_sets `TRUE` or `FALSE`. If `TRUE`, sets of SNPs for individual
##'   regions are returned, otherwise both sets and mapping information are
##'   returned. Default is `FALSE`.
##' @return A nested list containing following components:
##' - sets: a named list where each index represents a separate set of SNPs
##' - map: a data frame containing SNP mapping information
##' @export
map_snp_to_region <- function(info_snp, info_region,
                              extend_region_start = 20L,
                              extend_region_end = 20L,
                              only_sets=FALSE) {

  info_snp_name <- deparse(substitute(info_snp))

  is_df(info_snp)
  is_df(info_region)
  has_columns(info_snp, c("id", "chr", "pos"))
  has_columns(info_region, c("region.id", "chr", "start", "end"))
  is_nonnegative_number(extend_region_start, "extend_region_start")
  is_nonnegative_number(extend_region_end, "extend_region_end")

  if (anyDuplicated(info_snp$id) > 0L) {
    stop("SNP IDs in ", "'", info_snp_name, "'", " must be unique.",
         call. = FALSE)
  }

  if (anyDuplicated(info_region$region.id) > 0L) {
    warning(
      "Some region IDs have multiple intervals. ",
      "For those regions, SNPs mapped to any of these intervals ",
      "will be merged into a single set. ",
      "If you don't want this behavior, please assign unique IDs.",
      call. = FALSE
    )
  }

  snpdat <- data.table(
    id = as.character(info_snp$id),
    chr = as.character(info_snp$chr),
    start = as.integer(info_snp$pos),
    end = as.integer(info_snp$pos),
    key = c("chr", "start", "end")
  )

  regiondat <- data.table(
    region.id = as.character(info_region$region.id),
    chr = as.character(info_region$chr),
    start = as.integer(
      ifelse(info_region$start - extend_region_start * 1000 < 1L,
             1,
             info_region$start - extend_region_start * 1000)
    ),
    end = as.integer(info_region$end + extend_region_end * 1000),
    region.start = as.integer(info_region$start),
    region.end = as.integer(info_region$end),
    key = c("chr", "start", "end")
  )

  mapped <- foverlaps(snpdat, regiondat, type = "within")

  ## Unique function is necessary since non-overlapping regions with the same
  ## region.id can be overlapped by start:end position adjustment.
  ## e.g. RNU6-1: chr14 32202044:32202150 & 32203163:32203269
  ## (14:32183309) is doubly counted with region adj 20kb
  snp_sets <- lapply(split(mapped[!is.na(region.id), .(id, region.id)],
                           by = "region.id", keep.by = FALSE),
                     function(x) unique(unname(unlist(x))))

  if (only_sets) {
    list(sets = snp_sets)
  } else {
    setnames(mapped,
             old = c("start", "end", "i.start"),
             new = c("region.adj.start", "region.adj.end", "pos"))

    vkeep <- c("id", "chr", "pos",
               "region.id", "region.start", "region.end",
               "region.adj.start", "region.adj.end")

    list(sets = snp_sets,
         map = setDF(setorder(mapped[, vkeep, with = FALSE],
                              chr, region.start, region.end, pos,
                              na.last = TRUE)))
  }

 }
