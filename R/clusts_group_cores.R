#' @title clust_group_cores
#'
#' @description Takes the output of clust_groups, a data frame , and process
#' it to obtain the enrichment of genes in the core enrichment (or leading
#' edge analysis) within each cluster. The output is used in the functions
#' gm_enrichcores and gm_enrichreport.
#'
#' @param cg A data frame output from the GSEAmining cluster_groups function.
#' @param top An integer to choose the top most enriched genes to plot per
#' cluster. The default parameter are the top 3.
#'
#' @return
#' @export
#'
clust_group_cores <- function(cg,
                              top = 3) {
  # Obtain core enrichment terms of each gene set -----------------------------
  clust_lead <- cg %>%
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
    # eliminate tokens that are just numbers
    filter(!grepl('^\\d', lead_token)) %>%
    # group by cluster and count how many words
    group_by(Cluster, Enrichment) %>%
    count(lead_token, sort = TRUE) %>%
    # Select top n most common leading edge genes per cluster
    arrange(Cluster) %>%
    top_n(n = top)

}
