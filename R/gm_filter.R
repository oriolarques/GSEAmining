#' @title gm_filter: GSEAmining GSEA output filter
#'
#' @description Filters a data frame containing the results of GSEA analysis.
#'
#' @param df Data frame that contains at least three columns: an ID column for
#' the gene set names, a NES column with the normalized enrichment score and a
#' core_enrichment column containing the genes in the leading edge of each
#' gene set separated by '/'.
#' @param p.adj An integer to set the limit of the adjusted p-value (or false
#' discovery rate, FDR). Default value is 0.05
#' @param neg_NES A positive integer to set the limit of negative NES. Default
#' is 1.
#' @param pos_NES A positive integer to set the limit of positive NES. Default
#' is 1.
#'
#' @return A data frame.
#' @export
#'
#' @import dplyr
#'
#' @examples
#' data(genesets_sel)
#' gs.filt <- gm_filter(genesets_sel, p.adj = 0.05, neg_NES = 2.6, pos_NES = 2)
#'
gm_filter <- function(df,
                      p.adj = 0.05,
                      neg_NES = 1,
                      pos_NES = 1) {
  df %>%
    dplyr::arrange(desc(.data$NES)) %>%
    dplyr::filter((.data$p.adjust < p.adj & .data$NES < - neg_NES) |
                    (.data$p.adjust < p.adj & .data$NES > pos_NES)) %>%
    dplyr::select(.data$ID, .data$NES, .data$p.adjust, .data$core_enrichment)
}
