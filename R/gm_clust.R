#' @title gm_clust: GSEAmining cluster object
#' @description From the results of GSEA analysis, a data frame:
#' 1. Tokenizes the core-enrichment column.
#' 2. Cast the new table into a matrix.
#' 3. Returns and hclust object.
#'
#'
#' @param df Data frame that contains at least an ID column for the gene set
#' names and a core_enrichment column containing the genes in the leading edge
#' of each gene set separated by '/'.
#'
#'
#' @return A hclust class object with the clustering of the gene sets
#'  by the core enriched genes.
#'
#' @export
#'
#' @import dplyr
#' @import tidytext
#' @importFrom stats dist
#' @importFrom stats hclust
#'
#'
#'
#'
#'
gm_clust <- function(df) {

  stopifnot(is.data.frame(df))

  # Create a table where all leading edge genes are in rows -------------------
  gsea_tokens <- df %>%
    tidytext::unnest_tokens(output = lead_genes,
                            input = core_enrichment,
                            to_lower = FALSE) %>%
    group_by(ID) %>%
    count(lead_genes)

  # Casting tidytext data into a matrix ---------------------------------------
  cast_matrix <- cast_dtm(data = gsea_tokens,
                          document = ID,
                          term = lead_genes,
                          value = n)

  # Hierarchical clustering ---------------------------------------------------
  distance_matrix <- dist(cast_matrix, method = 'binary')

  hc <- hclust(distance_matrix)

}
