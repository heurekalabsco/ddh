test_that("get_gene_symbols_for_pathway returns gene symbols", {
  fake_get_data_object <- function (pathway_id, dataset_name) {
    id <- c('123')
    data_set <- c('universal_pathways')
    key <- c('gs_cat', 'gs_subcat', 'gene_symbol', 'gene_symbol', 'human_ensembl_gene')
    value <- c('C5', 'GO:BP', 'ROCK1', 'ROCK2', 'ENSG00000163347')
    tibble::as_tibble(data.frame(id, data_set, key, value))
  }
  local({
    mockr::local_mock(get_data_object = fake_get_data_object)
    gene_symbols <- get_gene_symbols_for_pathway('123')
    expect_equal(gene_symbols, c('ROCK1', "ROCK2"))
  })
})
