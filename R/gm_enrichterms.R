#' @title gm_enrichterms: GSEAmining enriched terms
#'
#' @description Takes the output of gm_clust, which is an hclust class object,
#' and plots gene set enriched terms as wordclouds. Two options are available,
#' either separate enrichments by clusters or plot them together in a single
#' plot.
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
#'
#' @return
#' @export
#'
#' @import stringr
#' @import ggwordcloud
#' @import ggplot2
#'
gm_enrichterms <- function(df,
                           hc,
                           clust = TRUE,
                           col_pos = 'red',
                           col_neg = 'blue') {

  # Get the cluster groups from  gm_clust object ------------------------------
  clust_groups <- data.frame(cluster=cutree(hc, h = 0.999)) %>%
    rownames_to_column('ID') %>%
    arrange(cluster) %>%
    left_join(df, by='ID')

  # Obtain the terms of each gene set -----------------------------------------
  stop_words <- readRDS(file = 'R/stop_words.rds')

  clust_groups_wordcloud <- clust_groups %>%
    mutate(enrichment = ifelse(test = NES > 0,
                               yes= 'pos',
                               no = 'neg')) %>%
    # separate the words of the geneset
    mutate(ID2 = str_replace_all(ID,
                                 pattern='_',
                                 replacement = ' ')) %>%
    # eliminate the first word of the geneset
    mutate(ID2 = word(.$ID2,
                      start = 2,
                      end = -1)) %>%
    # create word tokens (one row per word)
    tidytext::unnest_tokens(monogram,
                            ID2,
                            token = 'ngrams',
                            n = 1,
                            to_lower = FALSE) %>%
    # eliminate words like UP, DW, BY ...
    filter(!monogram %in% stop_words$word) %>%
    # eliminate words that are just numbers
    filter(!grepl('^\\d', monogram)) %>%
    # group by cluster and count how many words
    group_by(cluster, enrichment) %>%
    count(monogram, sort = TRUE)

  # Plot enriched terms wordclouds using ggwordcloud --------------------------
  plot <- ggplot(clust_groups_wordcloud,
                 aes(label = monogram,
                     size = n,
                     col = enrichment)) +
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
      facet_wrap(enrichment ~ cluster,
                 labeller = label_both)
  } else {
    plot +
      ggtitle('Gene Sets Enriched Terms')
  }

}
