##' Map SNPs to genes
##'
##' Annotate SNPs onto their neighboring genes (or arbitrary genomic regions) to
##' perform set-based association tests.
##'
##' @param info_snp A data frame with columns: "id", "chr", and "pos".
##' - id = SNP ID (e.g., rs numbers)
##' - chr = chromosome
##' - pos = base-pair position
##' @param info_gene A data frame with columns: "gene.id", "chr", "start",
##'   and "end".
##' - gene.id = gene ID (or identifier for genomic regions)
##' - chr = chromosome (must be the same chromosome coding scheme in `info_snp`)
##' - start = genomic start position
##' - end = genomic end position
##'
##' If a gene has multiple intervals, SNPs mapped to any of them will be merged
##' into a single set. Please assign unique IDs if you don't want this behavior.
##' @param extend_start A single non-negative integer, allowing for a certain kb
##'   window before the gene to be included. Default is 20 (= 20kb).
##' @param extend_end A single non-negative integer, allowing for a certain kb
##'   window after the gene to be included. Default is 20 (= 20kb).
##' @param only_sets If `TRUE`, only sets of SNPs for individual genes are
##'   returned. Otherwise, both sets and mapping information are returned.
##'   Default is `FALSE`.
##' @return A nested list containing following components:
##' - sets: a named list where each index represents a separate set of SNPs
##' - map: a data frame containing SNP mapping information
##' @examples
##' \dontshow{data.table::setDTthreads(1)}
##' ## GWAS summary statistics
##' head(exGWAS)
##'
##' ## Gene information data
##' head(gene.curated.GRCh37)
##'
##' ## Map SNPs to genes
##' snp_sets <- map_snp_to_gene(exGWAS, gene.curated.GRCh37)
##'
##' ## Better to use harmonized GWAS data for gene mapping
##' bfile <- system.file("extdata", "example.bed", package = "snpsettest")
##' x <- read_reference_bed(path = bfile)
##' hsumstats <- harmonize_sumstats(exGWAS, x)
##' snp_sets <- map_snp_to_gene(hsumstats, gene.curated.GRCh37)
##'
##' @export
map_snp_to_gene <- function(info_snp, info_gene,
                            extend_start = 20L,
                            extend_end = 20L,
                            only_sets=FALSE) {

  info_snp_name <- deparse(substitute(info_snp))

  is_df(info_snp)
  is_df(info_gene)
  has_columns(info_snp, c("id", "chr", "pos"))
  has_columns(info_gene, c("gene.id", "chr", "start", "end"))
  is_nonnegative_number(extend_start, "extend_start")
  is_nonnegative_number(extend_end, "extend_end")

  if (anyDuplicated(info_snp$id) > 0L) {
    stop("SNP IDs in ", "'", info_snp_name, "'", " must be unique.",
         call. = FALSE)
  }

  snpdat <- data.table(
    id = as.character(info_snp$id),
    chr = as.character(info_snp$chr),
    start = as.integer(info_snp$pos),
    end = as.integer(info_snp$pos),
    key = c("chr", "start", "end")
  )

  genedat <- data.table(
    gene.id = as.character(info_gene$gene.id),
    chr = as.character(info_gene$chr),
    start = as.integer(
      ifelse(info_gene$start - extend_start * 1000 < 1L,
             1,
             info_gene$start - extend_start * 1000)
    ),
    end = as.integer(info_gene$end + extend_end * 1000),
    gene.start = as.integer(info_gene$start),
    gene.end = as.integer(info_gene$end),
    key = c("chr", "start", "end")
  )

  mapped <- foverlaps(snpdat, genedat, type = "within")

  ## Unique function is necessary since non-overlapping genes with the same
  ## gene.id can be overlapped by start:end position adjustment.
  ## e.g. RNU6-1: chr14 32202044:32202150 & 32203163:32203269
  ## (14:32183309) is doubly counted with gene adj 20kb
  snp_sets <- lapply(split(mapped[!is.na(gene.id), .(id, gene.id)],
                           by = "gene.id", keep.by = FALSE),
                     function(x) unique(unname(unlist(x))))

  if (only_sets) {
    list(sets = snp_sets)
  } else {
    setnames(mapped,
             old = c("start", "end", "i.start"),
             new = c("gene.adj.start", "gene.adj.end", "pos"))

    vkeep <- c("id", "chr", "pos",
               "gene.id", "gene.start", "gene.end",
               "gene.adj.start", "gene.adj.end")

    list(sets = snp_sets,
         map = setDF(setorder(mapped[, vkeep, with = FALSE],
                              chr, gene.start, gene.end, pos,
                              na.last = TRUE)))
  }

 }
