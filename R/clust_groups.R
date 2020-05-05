#' @title clust_groups
#'
#' @description Takes the output of gm_clust, which is an hclust class object,
#' and returns a data frame that will be used in the rest of GSEAmining
#' functions gm_enrichreport, gm_enrichterms and gm_enrichcores.
#'
#' @param df Data frame that contains at least an ID column for the gene set
#' names and a core_enrichment column containing the genes in the leading edge
#' of each gene set separated by '/'.
#' @param hc The output of gm_clust, which is an hclust class object.
#'
#' @return A data.frame containing the cluster each gene set belongs to.
#' @export
#'
#'
clust_groups <- function(df,
                         hc) {

  # Get the cluster groups from  gm_clust object ------------------------------
  clust_groups <- data.frame(Cluster=cutree(hc, h = 0.999)) %>%
    rownames_to_column('ID') %>%
    arrange(.data$Cluster) %>%
    left_join(df, by='ID')

}
