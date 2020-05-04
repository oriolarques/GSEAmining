#' @title clust_group_terms
#'
#' @description Takes the output of clust_groups, a data frame , and process
#' it to obtain the enrichment of terms in gene sets names within each cluster.
#' The output is used in the functions gm_enrichterms and gm_enrichreport.
#'
#' @param cg A data frame output from the GSEAmining cluster_groups function.
#'
#' @return
#' @export
#'
clust_group_terms <- function(cg) {
  # Obtain the terms of each gene set -----------------------------------------
  stop.words <- stop_words()

  clust_groups_wordcloud <- cg %>%
    mutate(Enrichment = ifelse(test = NES > 0,
                               yes= 'Pos',
                               no = 'Neg')) %>%
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
    filter(!monogram %in% stop.words$word) %>%
    # eliminate words that are just numbers
    filter(!grepl('^\\d', monogram)) %>%
    # group by cluster and count how many words
    group_by(Cluster, Enrichment) %>%
    count(monogram, sort = TRUE)

}
