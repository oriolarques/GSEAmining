# GSEAmining
Make biological sense of Gene Set Enrichment Analysis outputs

Gene Set Enrichment Analysis 
([GSEA](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html)) 
is a very powerful and interesting computation method that allows an easy 
correlation between differential expressed genes to biological processes. 
Unfortunately, although it was designed to help researchers to interpret gene
expression data it can generate huge amounts of results whose biological 
meaning can be in turn difficult to interpret. 

Many available tools rely on the hierarchically structured Gene Ontology (GO) 
classification to reduce reundandcy in the results. However, due to the 
popularity of GSEA many more gene set collections, such as those in the 
Molecular Signatures Database 
([MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp)), 
are emerging. Since these collections 
are not organized as those in GO, their usage for GSEA do not always give a 
straightforward answer or, in other words, the results can be challenging to 
interpret with the currently available tools. For these reasons, GSEAmining 
was born to be an easy tool to create reproducible reports to help researchers 
make biological sense of GSEA outputs.

Given the results of GSEA, GSEAmining clusters the different gene sets based 
on the presence of the same genes in the leadind edge (or core) subset. 
Leading edge subset are those genes that contribute most to the enrichment 
score of each collection of genes or gene sets. For this reason, gene sets 
that participate in similar biological processes should share genes in common 
and in turn cluster together. After that, GSEAmining is able to identify and 
represent for each cluster two things: 

- The most enriched terms in the names of gene sets (as wordclouds)
- The most enriched genes in the leading edge subsets (as bar plots).

In each case, positive and negative enrichment of each gene sets are shown in 
different colors so it is easy to distinguish biological processes or genes
that may be of more interest in the particular case that is being studied.
