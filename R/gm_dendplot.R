#' @title gm_dendplot: GSEAmining dendrogram plotter
#'
#' @description Takes the output of gm_clust, which is an hclust class object,
#' and plots the dendrogram using the dendextend package.
#'
#' @param df Data frame that contains at least three columns: an ID column for
#' the gene set names, a NES column with the normalized enrichment score and a
#' core_enrichment column containing the genes in the leading edge of each
#' gene set separated by '/'.
#' @param hc The output of gm_clust, which is an hclust class object.
#' @param col_pos Color to represent the positively enriched gene sets. Default
#' is red.
#' @param col_neg Color to represent the negatively enriched gene sets. Default
#' is blue.
#' @param dend_len An integer that defines the length of the dendrogram. Default
#' value is 30. The closest to zero the longest the dendrogram.
#' @param rect A logical value indicating if rectangles should be drawn around
#' the clusters to help differentiating them. By default it is set to TRUE.
#' @param rect_len An integer to specify the length of the rectangle around the
#' cluster and the gene set label. Default is 2. The closest to zero the
#' smallest the rectangle.
#'
#' @return Invisibly returns a list with all the elements necessary to plot
#' a dendrogram.
#'
#' @export
#'
#' @import tibble
#' @import dendextend
#' @importFrom methods is
#'
#' @examples
#' data(genesets_sel)
#' gs.cl <- gm_clust(genesets_sel)
#' gm_dendplot(genesets_sel, gs.cl)
#'
gm_dendplot <- function(df,
                        hc,
                        col_pos = 'red',
                        col_neg ='blue',
                        dend_len = 30,
                        rect = TRUE,
                        rect_len = 2) {

  stopifnot(is.data.frame(df) | methods::is(hc, 'hclust'))

  # Convert hc in a dendrogram object to be visualized by dendextend
  dend <- hc %>% stats::as.dendrogram()

  # Color of the labels: Red=positively enriched / Blue = negatively enriched
  #   Get the order of the genesets in the cluster
  dend_order <- hc$labels[hc$order]

  # Create an object from geneset table and assign the color ------------------
  #   depending if NES < | > 0
  dend_label_col <- df %>%
    select(.data$ID, .data$NES) %>%
    mutate(color = ifelse(.data$NES > 0, col_pos, col_neg))

  #   Order the rows by the order of the genesets in the cluster
  dend_label_col <- dend_label_col[match(dend_order,
                                         dend_label_col$ID),]
  #   Keep only the value of color
  dend_label_col <- dend_label_col$color

  # Color of the clusters:
  #   Red=positively enriched / Blue = negatively enriched
  #   Identify to which cluster each geneset belongs to
  dend_clust_col <- data.frame(cluster=cutree(hc, h = 0.999))

  #   Get rownames to geneset column ------------------------------------------
  dend_clust_col <- dend_clust_col %>%
    rownames_to_column('geneset')

  # Join the cluster column to the geneset table ------------------------------
  dend_clust_col <- dend_clust_col%>%
    left_join(df, by=c('geneset'='ID')) %>%
    select(.data$geneset, .data$NES, .data$cluster) %>%
    # Define the color that each geneset has depending on NES
    mutate(color = ifelse(.data$NES > 0, col_pos, col_neg)) %>%
    select(.data$cluster, .data$color) %>%
    # Since positive and negative genesets cluster together eliminate duplicates
    unique() %>%
    # Sort the rows by cluster in descendant order
    # That's how dendextend reads clusters
    arrange(desc(.data$cluster))

  # Plot the clusters   -------------------------------------------------------
  graphics::par(mar = c(5.1,2,2,dend_len)) # set the margins
  dend %>%
    set("labels_cex", 0.45) %>%
    set('labels_col', dend_label_col) %>%
    set('branches_k_color', value = dend_clust_col$color,
        k = length(dend_clust_col[,1])) %>%
    plot(horiz = TRUE)

  # Add rectangles to differentiate better the clusters -----------------------
  if(rect == TRUE){
    dend %>%
      rect.dendrogram(k = length(dend_clust_col[,1]),
                      which = seq(from = 1,
                                  to = length((dend_clust_col[,1])),
                                  by = 2),
                      border = 0,
                      col = grDevices::rgb(red = 0.1,
                                           green =  0.2,
                                           blue =  0.4,
                                           alpha = 0.1),
                      lower_rect = -rect_len,
                      horiz = TRUE)
  }
}
