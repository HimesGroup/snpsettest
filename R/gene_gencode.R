##' Human gene information from the GENCODE GRCh37 version
##'
##' The gene information was extracted from the GENCODE release 19. This data
##' only included the following gene biotypes: protein-coding, Immunoglobulin
##' (Ig) variable chain and T-cell receptor (TcR) genes. Genes with
##' `NOVEL/Putative` status were excluded from the data.
##' @source <https://www.gencodegenes.org/human/release_19.html>
##' @format Data frame with columns
##' \describe{
##' \item{gene.id}{SNP ID.}
##' \item{chr}{chromosome.}
##' \item{start}{genomic start location (1-based).}
##' \item{end}{genomic end location.}
##' \item{strand}{genomic strand.}
##' \item{gene.name}{gene symbols mapped to the GENCODE genes.}
##' \item{gene.type}{gene biotypes in the GENCODE genes.}
##' }
##' @examples
##' ## Load data
##' data(gene.curated.GRCh37)
"gene.curated.GRCh37"

##' Human gene information from the GENCODE GRCh38 version
##'
##' The gene information was extracted from the GENCODE release 37. This data
##' only included the following gene biotypes: protein-coding, Immunoglobulin
##' (Ig) variable chain and T-cell receptor (TcR) genes.
##' @source <https://www.gencodegenes.org/human/release_37.html>
##' @format Data frame with columns
##' \describe{
##' \item{gene.id}{SNP ID.}
##' \item{chr}{chromosome.}
##' \item{start}{genomic start location (1-based).}
##' \item{end}{genomic end location.}
##' \item{strand}{genomic strand.}
##' \item{gene.name}{gene symbols mapped to the GENCODE genes.}
##' \item{gene.type}{gene biotypes in the GENCODE genes.}
##' }
##' @examples
##' ## Load data
##' data(gene.curated.GRCh38)
"gene.curated.GRCh38"
