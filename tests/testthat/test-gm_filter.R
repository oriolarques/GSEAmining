library("GSEAmining")
data(genesets_sel)

test_that("Test that gm_filter produces a data.frame with the 4 correct
          columns", {

  gs.filt <- gm_filter(genesets_sel,
                       p.adj = 0.05,
                       neg_NES = 2.6,
                       pos_NES = 2)

  expect_equal(class(gs.filt), 'data.frame')
  expect_equal(names(gs.filt), c("ID", "NES", "p.adjust", "core_enrichment"))

})
