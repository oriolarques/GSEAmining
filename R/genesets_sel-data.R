#' Selected gene sets as test
#'
#' Data that corresponds to GSEA analysis of differential expressed genes from
#' treated versus control samples in HGPalmer-PDX-P30 experiment. Differential
#' gene expression was obtained by using the oligo and limma R packages. GSEA
#' analysis was performed using the clusterProfiler R package using MSigDb
#' collections C2, C5 and Hallmarks.
#'
#' @docType data
#'
#' @usage data(genesets_sel)
#'
#' @format An object of class data.frame with 52 observations and 4 variables:
#' \describe{
#'     \item{ID}{Name of the gene set}
#'     \item{NES}{Normalized Enrichment Score}
#'     \item{p.adjust}{False discovery rate}
#'     \item{core_enrichment}{Genes that are in the leading edge subset}
#' }
#'
#' @references Arqués et al. Clinical Cancer Research. 2016 Feb 1;22(3):644-56.
#' doi: 10.1158/1078-0432.CCR-14-3081. Epub 2015 Jul 29.
#' (\href{https://clincancerres.aacrjournals.org/content/22/3/644.long}{Clinical Cancer Research})
#'
#' @source \href{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2446/}{ArrayExpress}
#'
#' @examples
#' data(genesets_sel)
#'
#' gs.cl <- gm_clust(genesets_sel)
#'
#' gm_dendplot(genesets_sel, gs.cl)
#'
#' gm_enrichterms(genesets_sel, gs.cl)
#'
#' gm_enrichcores(genesets_sel, gs.cl)
#'
#' \dontrun{gm_enrichreport(genesets_sel, gs.cl)}
#'
"genesets_sel"
