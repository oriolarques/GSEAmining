#' @title gm_enrichreport: GSEAmining Enrichment Report
#'
#' @description Takes the output of gm_clust, which is an hclust class object,
#' and creates a report in pdf that contains enriched terms and enriched core
#' genes in gene sets for each cluster. The results of each cluster are plotted
#' in an independent page.
#'
#' @param df Data frame that contains at least three columns: an ID column for
#' the gene set names, a NES column with the normalized enrichment score and a
#' core_enrichment column containing the genes in the leading edge of each
#' gene set separated by '/'.
#' @param hc The output of gm_clust, which is an hclust class object.
#' @param col_pos Color to represent positively enriched gene sets. Default
#' is red.
#' @param col_neg Color to represent negatively enriched gene sets. Default
#' is blue.
#' @param top An integer to choose the top most enriched genes to plot per
#' cluster. The default parameter are the top 3.
#' @param output A string to name the output pdf file.
#'
#' @return  Generates a pdf file.
#'
#' @export
#'
#' @importFrom gridExtra marrangeGrob
#'
#' @examples
#' #' data(genesets_sel)
#' gs.cl <- gm_clust(genesets_sel)
#' \dontrun{gm_enrichreport(genesets_sel, gs.cl)}
#'
gm_enrichreport <- function(df,
                            hc,
                            col_pos = 'red',
                            col_neg = 'blue',
                            top = 3,
                            output = 'gm_report') {

  stopifnot(is.data.frame(df) | !is(hc, 'hclust'))

  # 1. Get the cluster groups from  gm_clust object ---------------------------
  clust.groups <- clust_groups(df, hc)

  # 2. Obtain the terms of each gene set --------------------------------------
  clust.wc <- clust_group_terms(clust.groups)

  # 3. Obtain core enrichment terms of each gene set --------------------------
  clust.lead <- clust_group_cores(clust.groups, top)

  # 4. Put together gene set terms and core terms -----------------------------
  # Create a list containing ggplot wordcloud for enriched terms

  cgw <- list() # cgw: cluster.groups.wordclouds

  for(i in seq_along(unique(clust.groups$Cluster))){

    cgw[[i]] <- ggplot(clust.wc %>% filter(.data$Cluster == i),
                       aes(label = .data$monogram,
                           size = .data$n,
                           col = .data$Enrichment))+
      ggwordcloud::geom_text_wordcloud_area(eccentricity = 1,
                               area_corr_power = 1) +
      scale_size_area(max_size = 8)+ #max size the biggest word can have
      scale_color_manual(values = ifelse(
        clust.wc$Enrichment[clust.wc$Cluster==i][1] == 'Pos',
                                         col_pos, col_neg))+
      labs(title = paste0('Cluster ', clust.wc$Cluster[clust.wc$Cluster==i][1],
                    '\n Enrichment: ',
                    ifelse(clust.wc$Enrichment[clust.wc$Cluster==i][1] == 'Pos',
                      'Positive', 'Negative')),
           subtitle = 'Gene Set Terms Enriched') +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.07))
  }

  # Create a list containing ggplot barplots for core enrichment
  cgc <- list() # cgc: cluster.groups.core

  for(i in seq_along(unique(clust.groups$Cluster))){

    cgc[[i]] <- ggplot(clust.lead %>% filter(.data$Cluster == i),
                              aes(x = stats::reorder(.data$lead_token, .data$n),
                                  y = .data$n,
                                  fill = .data$Enrichment))+
      geom_bar(stat = 'identity') +
      ylim(0,max(clust.lead$n) + 1)+
      coord_flip()+
      scale_fill_manual(values = ifelse(clust.lead$Enrichment[
        clust.lead$Cluster==i][1] == 'Pos', col_pos, col_neg))+
      labs(title = '',
           subtitle = 'Top Genes in Leading Edge Analysis')+
      ylab('Num of Gene Appearances in Different Genesets')+
      theme_minimal() +
      theme(plot.subtitle = element_text(hjust = 0),
            axis.title.y = element_blank(),
            legend.position = 'none')
  }


  # Merge both tables (genesets wordcloud and leading edge of each geneset)
  # Intercalating each element wordcloud1, core1, wordcloud2, core2, ...
  plots <- c(rbind(cgw, cgc))

  # Save each combination of geneset wordcloud and leading edge in one page
  report <- marrangeGrob(grobs = plots,
                                   nrow = 2,
                                   ncol = 1,
                                   top = NULL)
  ggsave(paste0(output, '.pdf'), report)

}
