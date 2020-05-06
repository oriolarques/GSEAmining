#' @title clust_groups
#'
#' @description Takes the output of gm_clust, which is an hclust class object,
#' and returns a data frame that will be used in the rest of GSEAmining
#' functions gm_enrichreport, gm_enrichterms and gm_enrichcores.
#'
#' @param df Data frame that contains at least three columns: an ID column for
#' the gene set names, a NES column with the normalized enrichment score and a
#' core_enrichment column containing the genes in the leading edge of each
#' gene set separated by '/'.
#' @param hc The output of gm_clust, which is an hclust class object.
#'
#' @return A data.frame containing the cluster each gene set belongs to.
#'
#' @export
#'
#' @examples
#' data(genesets_sel)
#' gs.cl <- gm_clust(genesets_sel)
#' clust.groups <- clust_groups(genesets_sel, gs.cl)
#'
clust_groups <- function(df,
                         hc) {

  # Get the cluster groups from  gm_clust object ------------------------------
  clust_groups <- data.frame(Cluster=cutree(hc, h = 0.999)) %>%
    rownames_to_column('ID') %>%
    arrange(.data$Cluster) %>%
    left_join(df, by='ID')

}
