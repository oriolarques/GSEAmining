gm_enrichall <- function(df,
                         hc) {
  #### Put together geneset wordcloud and leading edge analysis ---------------------------------
  # Create a list containing ggplot wordcloud for genesets
  mining_pgs <- list()
  for(i in 1:length(unique(mining_clusters$cluster))){

    mining_pgs[[i]] <- ggplot(mining_clusters_wordcloud %>% filter(cluster == i),
                              aes(label=monogram, size=n, col = enrichment))+
      geom_text_wordcloud_area(eccentricity = 1, area_corr_power = 1) +
      scale_size_area(max_size = 8)+ #max size the biggest word can have
      scale_color_manual(values = ifelse(mining_clusters_wordcloud$enrichment[mining_clusters_wordcloud$cluster==i][1] == 'pos',
                                         'red', 'blue'))+
      labs(title = paste0('Cluster ', mining_clusters_wordcloud$cluster[mining_clusters_wordcloud$cluster==i][1],
                          '\n Enrichment: ', ifelse(mining_clusters_wordcloud$enrichment[mining_clusters_wordcloud$cluster==i][1] == 'pos', 'Positive', 'Negative')),
           subtitle = 'Gene Set Terms Enriched') +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.07))
  }

  # Create a list containing ggplot barplots for leading edge
  mining_ple <- list()
  for(i in 1:length(unique(mining_clusters$cluster))){

    mining_ple[[i]] <- ggplot(mining_clusters_leading_wordcloud %>% filter(cluster == i),
                              aes(x=reorder(lead_token, n), y=n, fill = enrichment))+
      geom_bar(stat = 'identity') +
##############!!!!!!!!!!!!!!!!!!!ylim(0,42)+
      coord_flip()+
      scale_fill_manual(values = ifelse(mining_clusters_leading_wordcloud$enrichment[mining_clusters_leading_wordcloud$cluster==i][1] == 'pos',
                                        'red', 'blue'))+
      labs(title = '', subtitle = 'Top Genes in Leading Edge Analysis')+
      ylab('Num of Gene Appearances in Different Genesets')+
      theme_minimal() +
      theme(plot.subtitle = element_text(hjust = 0),
            axis.title.y = element_blank(),
            legend.position = 'none')
  }


  # Merge both tables (genesets wordcloud and leading edge of each geneset)
  # Intercalating each element wordcloud1, leadedge1, wordcloud2, leadedge2, ...
  mining_p <- c(rbind(mining_pgs, mining_ple))


  # Plot side by side each geneset wordcloud with each leading edge
  #grid.arrange(grobs=mining_p, nrow=3)


  # Save each combination of geneset wordcloud and leading edge in one page
  mining_finalplot <- marrangeGrob(grobs = mining_p, nrow = 2, ncol = 1, top = NULL)
  ggsave('finalplot2.pdf', mining_finalplot)

}
