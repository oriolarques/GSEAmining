#' @title gm_clust: GSEAmining cluster object
#'
#' @description Takes the output of gm_filter or a data frame that with the
#' results of GSEA analysis and returns and hclust object that can be plotted
#' using the gm_dendplot function.
#'
#'
#' @param df Data frame that contains at least three columns: an ID column for
#' the gene set names, a NES column with the normalized enrichment score and a
#' core_enrichment column containing the genes in the leading edge of each
#' gene set separated by '/'.
#'
#'
#' @return An object of class hclust that contains the clustering of the gene
#' sets by the core enriched genes.First a distance matrix is calculated
#' using the 'binary' method and then a cluster with the 'complete' method is
#' created.
#'
#' @export
#'
#' @import dplyr
#' @import tidytext
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom rlang .data
#'
#' @examples
#' data(genesets_sel)
#' gs.cl <- gm_clust(genesets_sel)
#'
#'
gm_clust <- function(df) {

  stopifnot(is.data.frame(df))

  # Create a table where all leading edge genes are in rows -------------------
  gsea_tokens <- df %>%
    tidytext::unnest_tokens(output = .data$lead_genes,
                            input = .data$core_enrichment,
                            to_lower = FALSE) %>%
    group_by(.data$ID) %>%
    count(.data$lead_genes)

  # Casting tidytext data into a matrix ---------------------------------------
  cast_matrix <- cast_dtm(data = gsea_tokens,
                          document = .data$ID,
                          term = .data$lead_genes,
                          value = n)

  # Hierarchical clustering ---------------------------------------------------
  distance_matrix <- dist(cast_matrix, method = 'binary')

  hc <- hclust(distance_matrix)

}
