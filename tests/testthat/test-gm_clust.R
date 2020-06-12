library("GSEAmining")
data(genesets_sel)

test_that("Test that gm_clust produces an hclust object", {

  gs.filt <- gm_filter(genesets_sel,
                       p.adj = 0.05,
                       neg_NES = 2.6,
                       pos_NES = 2)
  gs.cl <- gm_clust(gs.filt)

  expect_equal(class(gs.cl), 'hclust')
  expect_equal(length(gs.cl), 7)

})
