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
#' @param neg_NES An integer to set the limit of negative NES. Default is NULL.
#' @param pos_NES An integer to set the limit of positive NES. Default is NULL.
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
                      neg_NES = NULL,
                      pos_NES = NULL) {
  df %>%
    dplyr::arrange(desc(NES)) %>%
    dplyr::filter((p.adjust < p.adj & NES < - neg_NES) |
                    (p.adjust < p.adj & NES > pos_NES)) %>%
    dplyr::select(ID, NES, p.adjust, core_enrichment)
}
