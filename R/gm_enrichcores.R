#' @title gm_enrichcores: GSEAmining core enrichment genes
#'
#' @description Takes the output of gm_clust, which is an hclust class object,
#' and plots the top n genes in core enrichment (leading edge analysis).
#' Two options are available, either separate barplots by clusters or all
#' together in one plot.
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
#' @param top An integer to choose the top most enriched genes to plot per
#' cluster. The default parameter are the top 3.
#'
#' @return Returns a ggplot object.
#'
#' @export
#'
#' @importFrom grDevices rgb
#' @importFrom graphics par
#'
#' @examples
#' data(genesets_sel)
#' gs.cl <- gm_clust(genesets_sel)
#' gm_enrichcores(genesets_sel, gs.cl)
#'
gm_enrichcores <- function(df,
                           hc,
                           clust = TRUE,
                           col_pos = 'red',
                           col_neg = 'blue',
                           top = 3) {

  stopifnot(is.data.frame(df) | methods::is(hc, 'hclust'))

  # Get the cluster groups from  gm_clust object ------------------------------
  clust.groups <- clust_groups(df, hc)

  # Obtain core enrichment terms of each gene set -----------------------------
  clust.lead <- clust_group_cores(clust.groups, top)

  # Bar plot of top n core enrichment gens ------------------------------------
  plot <- ggplot(clust.lead,
                 aes(x = stats::reorder(.data$lead_token, .data$n),
                     y = .data$n,
                     fill = .data$Enrichment)) +
    geom_bar(stat = 'identity') +
    coord_flip() +
    scale_fill_manual(values = c(col_neg, col_pos)) +
    scale_y_continuous(expand = c(0, 0)) +
    ylab('Num of Gene Appearances in Different Genesets') +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          panel.grid.major.y = element_blank(),
          strip.text = element_text(face = 'bold'))

  if (clust == TRUE) {
    plot +
      ggtitle('Top Genes in Leading Edge Analysis by Cluster') +
      theme(legend.position = 'none') +
      facet_wrap(Enrichment ~ Cluster,
                 labeller= label_both,
                 scales = 'free_y')
  } else {
    plot +
      ggtitle('Top Genes in Leading Edge Analysis')
  }

}
