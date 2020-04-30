#' @title gm_enrichcores: GSEAmining core enrichment genes
#'
#' @description Takes the output of gm_clust, which is an hclust class object,
#' and plots the top n genes in core enrichment (leading edge analysis).
#' Two options are available, either separate barplots by clusters or all
#' together in one plot.
#'
#' @param df Data frame that contains at least an ID column for the gene set
#' names and a core_enrichment column containing the genes in the leading edge
#' of each gene set separated by '/'.
#' @param hc The output of gm_clust, which is an hclust class object.
#' @param clust A logical value indicating if wordclouds should be separated by
#' clusters or not. Default value is TRUE.
#' @param col_pos Color to represent the positively enriched gene sets. Default
#' is red.
#' @param col_neg Color to represent the negatively enriched gene sets. Default
#' is blue.
#' @param top An integer to choose he top most enriched genes to plot per
#' cluster. The default parameter are the top 3.
#'
#' @return
#' @export
#'
#'
#'
gm_enrichcores <- function(df,
                           hc,
                           clust = TRUE,
                           col_pos = 'red',
                           col_neg = 'blue',
                           top = 3) {

  stopifnot(is.data.frame(df) | class(hc) != 'hclust')

  # Get the cluster groups from  gm_clust object ------------------------------
  clust_groups <- data.frame(Cluster=cutree(hc, h = 0.999)) %>%
    rownames_to_column('ID') %>%
    arrange(Cluster) %>%
    left_join(df, by='ID')

  # Obtain core enrichment terms of each gene set -----------------------------
  clust_lead <- clust_groups %>%
    mutate(Enrichment = ifelse(NES > 0, 'Pos', 'Neg')) %>%
    # separate the words of the core enrichment
    mutate(core_enrichment = str_replace_all(core_enrichment,
                                             pattern='/',
                                             replacement = ' ')) %>%
    # create gene word tokens (one row per word)
    tidytext::unnest_tokens(lead_token,
                            core_enrichment,
                            token = 'ngrams',
                            n = 1,
                            to_lower = FALSE) %>%
    # group by cluster and count how many words
    group_by(Cluster, Enrichment) %>%
    count(lead_token, sort = TRUE) %>%
    # Select top n most common leading edge genes per cluster
    arrange(Cluster) %>%
    top_n(n = top)


  # Bar plot of top n core enrichment gens ------------------------------------
  plot <- ggplot(clust_lead,
                 aes(x = reorder(lead_token, n),
                     y = n,
                     fill = Enrichment)) +
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
