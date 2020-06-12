library("GSEAmining")
data(genesets_sel)

test_that("Test that clust_groups, clust_group_terms and clust_group_cores
          produce data.frames with necessary variables", {

  gs.filt <- gm_filter(genesets_sel,
                       p.adj = 0.05,
                       neg_NES = 2.6,
                       pos_NES = 2)
  gs.cl <- gm_clust(gs.filt)

  cg <- clust_groups(genesets_sel, gs.cl)
  cg_t <- clust_group_terms(cg)
  cg_c <- clust_group_cores(cg)


  expect_equal(class(cg), 'data.frame')
  expect_equal(names(cg),
               c("ID", "Cluster", "NES", "p.adjust", "core_enrichment"))
  expect_equal(names(cg_t),
               c("Cluster", "Enrichment", "monogram", "n"))
  expect_equal(names(cg_c),
               c("Cluster", "Enrichment", "lead_token", "n"))

})
