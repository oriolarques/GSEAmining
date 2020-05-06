#' @title gm_enrichterms: GSEAmining enriched terms
#'
#' @description Takes the output of gm_clust, which is an hclust class object,
#' and plots gene set enriched terms as wordclouds. Two options are available,
#' either separate enrichments by clusters or plot them together in a single
#' plot.
#'
#' @param df Data frame that contains at least three columns: an ID column for
#' the gene set names, a NES column with the normalized enrichment score and a
#' core_enrichment column containing the genes in the leading edge of each
#' gene set separated by '/'.
#' @param hc The output of gm_clust, which is an hclust class object.
#' @param clust A logical value indicating if wordclouds should be separated by
#' clusters or not. Default value is TRUE.
#' @param col_pos Color to represent positively enriched gene sets. Default
#' is red.
#' @param col_neg Color to represent negatively enriched gene sets. Default
#' is blue.
#'
#' @return Returns a ggplot object.
#'
#' @export
#'
#' @import stringr
#' @import ggwordcloud
#' @import ggplot2
#'
#' @examples
#' data(genesets_sel)
#' gs.cl <- gm_clust(genesets_sel)
#' gm_enrichterms(genesets_sel, gs.cl)
#'
gm_enrichterms <- function(df,
                           hc,
                           clust = TRUE,
                           col_pos = 'red',
                           col_neg = 'blue') {

  stopifnot(is.data.frame(df) | methods::is(hc, 'hclust'))

  # Get the cluster groups from  gm_clust object ------------------------------
  clust.groups <- clust_groups(df, hc)

  # Obtain the terms of each gene set -----------------------------------------
  clust.groups.wordcloud <- clust_group_terms(clust.groups)

  # Plot enriched terms wordclouds using ggwordcloud --------------------------
  plot <- ggplot(clust.groups.wordcloud,
                 aes(label = .data$monogram,
                     size = .data$n,
                     col = .data$Enrichment)) +
    # rm_outside option eliminates words that do not fit in the image
    ggwordcloud::geom_text_wordcloud_area(eccentricity = 1,
                                          area_corr_power = 1,
                                          rm_outside = TRUE) +
    scale_color_manual(values = c(col_neg, col_pos)) +
    scale_size_area(max_size = 8) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          strip.text = element_text(face = 'bold'))

  # Plot text title and facets depending on clustering parameter --------------
  if (clust == TRUE) {
    plot +
      ggtitle('Gene Sets Enriched Terms  by Cluster') +
      facet_wrap(Enrichment ~ Cluster,
                 labeller = label_both)
  } else {
    plot +
      ggtitle('Gene Sets Enriched Terms')
  }

}
