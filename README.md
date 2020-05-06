# GSEAmining
Make biological sense of Gene Set Enrichment Analysis outputs

Gene Set Enrichment Analysis (GSEA) is a very powerful and interesting 
computation method that allows and easy correlation between differential 
expressed genes to biological processes.Unfortunately, although it was 
designed to help researchers to interpret gene expression data it can 
generate huge amounts of results whose biological meaning can be in turn 
difficult to interpret. 
Many available tools rely on hierarchically structured Gene Ontology (GO) 
classification to reduce reundandcy in the results. However, due to the 
popularity of GSEA many more gene set collections, such as those in the 
Molecular Signatures Database (MSigDB), are emerging. Since these collections 
are not organized as those in GO, their usage for GSEA do not always give a 
straightforward answer or, in other words, the results can be challenging to 
interpret with the currently available tools. For these reasons, GSEAmining 
was born to be an easy tool to create reproducible reports to help researchers
make biological sense of GSEA outputs.
