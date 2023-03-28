#' Stop words.
#' Eliminates the first word of the gene sets from MSigDb that relate to the
#' origin of the gene set. Additionally it eliminates words that do not add
#' a lot of significance such as prepositions or adverbs among others.
#'
#' @return Returns a tibble with 2 variables.
#'
#' @export
#'
#'
stop_words <- function() {
  return(tribble(
  ~word, ~lexicon,
  'REACTOME', "CUSTOM",
  'HALLMARK', "CUSTOM",
  'KEGG', "CUSTOM",
  'BIOCARTA', "CUSTOM",
  'PID', "CUSTOM",
  'TARGETS', "CUSTOM",
  'RESPONSE', "CUSTOM",
  'PATHWAY', "CUSTOM",
  'GO', "CUSTOM",
  'GOBP', "CUSTOM",
  'GOCC', "CUSTOM",
  'GOMF', "CUSTOM",
  'OF', "CUSTOM",
  'BY', "CUSTOM",
  'AND', "CUSTOM",
  'THE', "CUSTOM",
  'IN', "CUSTOM",
  'TO', "CUSTOM",
  'UP', "CUSTOM",
  'DN', "CUSTOM",
  'OF', "CUSTOM",
  'VS', "CUSTOM",
  'OR', "CUSTOM",
  'AT', "CUSTOM",
  'REGULATION', "CUSTOM",
  'PROCESS', "CUSTOM",
  'POSITIVE', "CUSTOM",
  'NEGATIVE', "CUSTOM",
  'FROM', "CUSTOM",
  'INTO', "CUSTOM",
  'UNKNOWN', "CUSTOM",
  'MODULE', "CUSTOM",
  '1HR', "CUSTOM",
  '2HR', "CUSTOM",
  '4HR', "CUSTOM",
  '6HR', "CUSTOM",
  '8HR', "CUSTOM",
  '10HR', "CUSTOM",
  '12HR', "CUSTOM",
  '14HR', "CUSTOM",
  '16HR', "CUSTOM",
  '18HR', "CUSTOM",
  '20HR', "CUSTOM",
  '24HR', "CUSTOM",
  '48HR', "CUSTOM",
  '96HR', "CUSTOM"))
}
