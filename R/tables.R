
# PATHWAY TABLES -----
#' Pathway List Table
#'
#' \code{make_pathway_list} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a sub-table of gene sets and pathways that contain your gene query
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a pathway list table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_pathway_list(input = list(type = 'gene', content = 'ROCK1'))
#' make_pathway_list(input = list(type = 'gene', content = c('ROCK1', 'ROCK2')))
#' \dontrun{
#' make_pathway_list(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_pathway_list <- function(input = list()) {
  make_pathway_list_raw <- function() {
    gene_pathways <-
      get_data_object(object_name = input$content,
                      dataset_name = "universal_gene_pathways",
                      pivotwider = TRUE) %>%
      dplyr::select(-any_of(c("id", "data_set")))
    return(gene_pathways)
  }
  #error handling
  tryCatch(make_pathway_list_raw(),
           error = function(e){
             message(e)
           })
}

#' Pathway Genes Table
#'
#' \code{make_pathway_genes} returns a table of genes in a queried pathway
#'
#' This is a table function that takes a gene_set id in the query slot and returns a sub-table of genes that your gene set contains
#'#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_pathway_genes(input = list(type = 'gene', subtype = 'pathway', query = '16769'))
make_pathway_genes <- function(input = list()){
  make_pathway_genes_raw <- function() {
    gene_pathways <-
      get_data_object(object_name = input$query,
                      dataset_name = "universal_pathways",
                      pivotwider = TRUE) %>%
      dplyr::select(gene_symbol) %>%
      tidyr::unnest(cols = "gene_symbol")
    gene_pathway_table <-
      gene_pathways %>%
      # dplyr::filter(gene_symbol %in% c("ROCK1", "ROCK2")) %>% #for testing
      dplyr::mutate(gene_description = purrr::map_chr(gene_symbol, ~ make_summary_gene(input = list(content = .), var = "approved_name")))

    return(gene_pathway_table)
  }
  #error handling
  tryCatch(make_pathway_genes_raw(),
           error = function(e){
             message(e)
           })
}

#' Compound Table
#'
#' \code{make_compound_table} returns an image of ...
#'
#' This is a table function that takes a compound name and returns a table. It is commented out until we revisit compouter as a query feather.
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a compound Table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_compound_table(input = list(content = "aspirin"), top = TRUE)
#' \dontrun{
#' make_compound_table(input = list(content = "aspirin"), top = FALSE)
#' }
make_compound_table <- function(input = list(),
                                # data_compound_prism_cor_nest = compound_prism_cor_nest,
                                # data_compound_prism_names = compound_prism_names,
                                # data_universal_stats_summary = universal_stats_summary,
                                top = TRUE) {
  NULL
  # make_compound_table_raw <- function() {
  #   prism_cor_upper <- get_stats(data_set = "prism", var = "upper")
  #   prism_cor_lower <- get_stats(data_set = "prism", var = "lower")
  #
  #   table_complete <-
  #     data_compound_prism_cor_nest %>%
  #     dplyr::filter_all(dplyr::any_vars(fav_drug %in% input$content)) %>%
  #     tidyr::unnest("data") %>%
  #     dplyr::ungroup() %>%
  #     {if (top == TRUE) dplyr::filter(., r2 > prism_cor_upper) else dplyr::filter(., r2 < prism_cor_lower)} %>%
  #     dplyr::arrange(dplyr::desc(r2)) %>%
  #     dplyr::left_join(data_compound_prism_names, by = "name") %>%
  #     dplyr::select(1:4)
  #   return(table_complete)
  # }
  # #error handling
  # tryCatch(make_compound_table_raw(),
  #          error = function(e){
  #            #make empty table equivalent to returned table
  #            return(data_prism_cor_nest %>%
  #                     tidyr::unnest("data") %>%
  #                     dplyr::ungroup() %>%
  #                     dplyr::left_join(data_compound_prism_names, by = "name") %>%
  #                     dplyr::select(1:4) %>%
  #                     dplyr::slice(0))
  #          }
  # )
}

# PUBMED TABLE -----
#' Pubmed Table
#'
#' \code{make_pubmed_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a pubmed Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a pubmed table..
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_pubmed_table(input = list(type = 'gene', content = 'ROCK1'))
#' \dontrun{
#' make_pubmed_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_pubmed_table <- function(input = list()) {
  make_pubmed_table_raw <- function() {
    pubmed_table <-
      get_data_object(object_name = input$content,
                      dataset_name = "universal_pubmed",
                      pivotwider = TRUE) %>%
      dplyr::mutate(year = as.numeric(year),
                    pmid = as.numeric(pmid)) %>%
      dplyr::arrange(as.numeric(pmid)) %>%
      dplyr::select(id, pmid, year, pmcid)
    return(pubmed_table)
  }
  #error handling
  tryCatch(make_pubmed_table_raw(),
           error = function(e){
             message(e)
           })
}

# CELL ANATOGRAM TABLES -----
#' Cellanatogram Table
#'
#' \code{make_cellanatogram_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a cellanatogram Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a cellanatogram Table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_cellanatogram_table(input = list(type = 'gene', content = 'ROCK2'))
#' make_cellanatogram_table(input = list(type = 'gene', content = c('ROCK1', 'ROCK2')))

#' \dontrun{
#' make_cellanatogram_table(input = list(type = 'gene', content = 'ROCK2'))
#' }
make_cellanatogram_table <- function(input = list()) {
  make_cellanatogram_table_raw <- function() {

    cellanatogram_table <-
      get_data_object(object_names = input$content,
                      dataset_name = "gene_subcell",
                      pivotwider = TRUE) %>%
      dplyr::mutate(value = round(as.numeric(value), 1)) %>%
      dplyr::arrange(dplyr::desc(value))
    return(cellanatogram_table)
  }
  #error handling
  tryCatch(make_cellanatogram_table_raw(),
           error = function(e){
             message(e)
           })
}

# EXPRESSION TABLES -----
#' Expression Table
#'
#' \code{make_expression_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a expression Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a expression Table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_expression_table(input = list(type = 'gene', content = 'ROCK1'))
#' make_expression_table(input = list(type = 'gene', content = 'ROCK1'), var = "protein")
#' make_expression_table(input = list(type = 'gene', content = c('ROCK1', 'ROCK2')))
#' \dontrun{
#' make_expression_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_expression_table <- function(input = list(),
                                  var = "gene") { #you are so slow
  make_expression_table_raw <- function() {
    data_universal_expression_long <-
      get_data_object(object_names = input$content,
                      dataset_name = "universal_expression_long",
                      pivotwider = TRUE) %>%
      dplyr::mutate(across(contains(c("expression")), as.numeric)) %>%
      dplyr::select(-data_set)

    cell_expression_names <- get_content("cell_expression_names", dataset = TRUE)

    if (var == "gene") {
      table_data <-
        data_universal_expression_long %>%
        dplyr::select(-"protein_expression") %>%
        dplyr::rename("expression_var" = "gene_expression")
    } else if (var == "protein") {
      table_data <-
        data_universal_expression_long %>%
        dplyr::select(-"gene_expression") %>%
        dplyr::rename("expression_var" = "protein_expression")
    } else {
      stop("delcare your variable")
    }
    if (input$type == "gene") {
      expression_table <-
        table_data %>%
        dplyr::arrange(dplyr::desc(.[3])) %>%
        tidyr::pivot_wider(names_from = id, values_from = expression_var) %>%
        dplyr::left_join(cell_expression_names, by = "depmap_id") %>%
        dplyr::select(-depmap_id) %>%
        dplyr::select(cell_line, lineage, lineage_subtype, dplyr::everything()) %>%
        dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
        dplyr::rename("Cell Line" = "cell_line", "Lineage" = "lineage", "Subtype" = "lineage_subtype")
    } # else if (input$type == "cell") {
    # table_data %>%
    #   dplyr::arrange(dplyr::desc(.[3])) %>%
    #   dplyr::left_join(data_cell_expression_names, by = "X1") %>%
    #   dplyr::select(-X1, -lineage, -lineage_subtype) %>%
    #   dplyr::filter_all(dplyr::any_vars(cell_line %in% input$content)) %>%
    #   tidyr::pivot_wider(names_from = cell_line, values_from = expression_var) %>%
    #   dplyr::select(gene, dplyr::everything()) %>%
    #   dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
    #   dplyr::rename("Gene" = "gene") %>%
    #   dplyr::left_join(data_universal_gene_summary %>%
    #                      dplyr::select(Gene = approved_symbol,
    #                                    `Gene Name` = approved_name),
    #                    by = "Gene") %>%
    #   dplyr::relocate(`Gene Name`, .after = Gene)
    # }
    return(expression_table)
  }
  #error handling
  tryCatch(make_expression_table_raw(),
           error = function(e){
             message(e)
           })
}

# HUMAN ANATOGRAM TABLES -----
#' Human Anatogram Table
#'
#' \code{make_humananatogram_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a humananatogram Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a human anatogram Table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_humananatogram_table(input = list(type = 'gene', content = 'ROCK1'))
#' make_humananatogram_table(input = list(type = 'gene', content = c('ROCK1', 'ROCK2')))
#' \dontrun{
#' make_humananatogram_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_humananatogram_table <- function(input = list()) {
  make_humananatogram_table_raw <- function() {
    humananatogram_table <-
      get_data_object(object_names = input$content,
                      dataset_name = "gene_tissue",
                      pivotwider = TRUE) %>%
      dplyr::mutate(across(contains(c("value")), as.numeric)) %>%
      dplyr::select(id, organ, value) %>%
      dplyr::filter(!is.na(value)) %>%
      dplyr::mutate(organ = stringr::str_replace_all(organ, "_", " "),
                    organ = stringr::str_to_title(organ)) %>%
      dplyr::arrange(dplyr::desc(value)) %>%
      tidyr::pivot_wider(names_from = id, values_from = value, names_sort = TRUE) %>%
      dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
      dplyr::rename("Organ" = "organ")
    return(humananatogram_table)
  }
  #error handling
  tryCatch(make_humananatogram_table_raw(),
           error = function(e){
             message(e)
           })
}
## GET CLUSTER  -----------------------------------------------
#' Clustering Table
#'
#' This is a helper function that takes a gene name and return a vector containing the cluster numbers for input
#'
#' @param input Expecting a list containing a content variable.
#' @return Returns a vector containing the cluster number.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' get_cluster(input = list(type = 'gene', content = 'ROCK1'))
#' get_cluster(input = list(type = 'gene', content = c('ROCK1', 'ROCK2')))
#' \dontrun{
#' get_cluster(input = list(type = 'gene', content = 'ROCK1'))
#' }
get_cluster <- function(input = list()){
  cluster <- get_data_object(object_names = input$content,
                             dataset_name = "gene_signature_clusters") %>%
    dplyr::filter(key == "clust") %>%
    dplyr::pull(value) %>%
    unique()
  return(cluster)
}

## PROTEIN CLUSTER TABLE -----------------------------------------------
#' Clustering Table
#'
#' This is a table function that takes a gene name and returns a clustering Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a clustering Table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_signature_clusters_table(input = list(type = 'gene', content = 'ROCK1'))
#' \dontrun{
#' make_signature_clusters_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_signature_clusters_table <- function(input = list(),
                                          ...) {
  make_signature_clusters_table_raw <- function() {

    gene_signature_clusters <- dplyr::as_tibble(get_content("gene_signature_clusters", dataset = TRUE)) # allows all genes, coordinates, and cluster

    #vec of clusters in query
    query_clust <-
      ddh::get_cluster(input)

    cluster_table <-
      gene_signature_clusters %>%
      dplyr::filter(clust %in% query_clust) %>%
      dplyr::select(Gene = id, Name = protein_name, Cluster = clust) %>%
      dplyr::as_tibble()

    return(cluster_table)
  }

  #error handling
  tryCatch(make_signature_clusters_table_raw(),
           error = function(e){
             message(e)
           })
}

## PROTEIN CLUSTER ENRICHMENT TABLE -----------------------------------------------
#' Clustering Enrichment Table
#'
#' \code{make_clustering_enrichment_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a clustering enrichment Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a clustering enrichment Table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_clustering_enrichment_table(input = list(type = 'gene', content = 'ROCK1'))
#' make_clustering_enrichment_table(input = list(type = 'gene', content = 'ROCK1'), ontology = "MF")
#' make_clustering_enrichment_table(input = list(type = 'gene', content = c('ROCK1', 'ROCK2')))
#' \dontrun{
#' make_clustering_enrichment_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_clustering_enrichment_table <- function(input = list(),
                                             ontology = "BP",
                                             filter_noise = FALSE) {

  make_clustering_enrichment_table_raw <- function() {

    gene_signature_clusters <- ddh::get_content("gene_signature_clusters", dataset = TRUE) #allows all genes, coordinates, and cluster
    gene_signature_cluster_enrichment <- ddh::get_content("gene_signature_cluster_enrichment", dataset = TRUE)

    if(filter_noise) {
      gene_signature_clusters <-
        gene_signature_clusters %>%
        dplyr::filter(clust != 0)
    }

    query_clust <- ddh::get_cluster(input = input)

    if(length(query_clust) != 1) {
      stop("Select only one cluster")
    }

    enrichment_table <-
      gene_signature_cluster_enrichment %>%
      dplyr::filter(clust %in% query_clust) %>%
      dplyr::filter(ont %in% ontology) %>%
      dplyr::select(-all_of(c("ont", "clust")))

    if(nrow(enrichment_table) == 0) {
      stop("No enriched terms for this cluster and ontology...")
    }
    return(enrichment_table)
  }

  #error handling
  tryCatch(make_clustering_enrichment_table_raw(),
           error = function(e){
             message(e)
           })
}

## 3D STRUCTURE TABLE --------------------------------------------------------------------
#' Plot for 3D protein structure table
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_structure3d_table(input = list(content = 'ROCK1'))
#' #' \dontrun{
#' make_structure3d_table(input = list(content = 'ROCK1'))
#' }
make_structure3d_table <- function(data_gene_uniprot_pdb_table = gene_uniprot_pdb_table,
                                   data_universal_proteins = universal_proteins,
                                   input = list()) {
  make_structure3d_table_raw <- function() {
    table_complete <-
      ddh::get_data_object(object_name = input$content,
                           dataset_name = "gene_pdb_table",
                           pivotwider = TRUE) %>%
      dplyr::select(tidyselect::any_of(c("id", "pdb", "title", "doi", "organism", "expression_system")))

    return(table_complete)
  }

  #error handling
  tryCatch(make_structure3d_table_raw(),
           error = function(e){
             return()
           })
}

# DEPENDENCY TABLES -----
#' Dependency Table
#'
#' \code{make_dep_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a dep Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a dep Table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_dep_table(input = list(type = 'gene', content = 'ROCK1'))
#' make_dep_table(input = list(type = 'gene', content = c('ROCK1', 'ROCK2')))
#' make_dep_table(input = list(type = 'compound', content = 'aspirin'))
#' make_dep_table(input = list(type = 'cell', content = 'HEPG2'))
#' \dontrun{
#' make_dep_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_dep_table <- function(input = list()#,
                           # leave these here until you go back and fix cell and drug
                           # data_universal_prism_long = universal_prism_long,
                           # data_compound_prism_names = compound_prism_names,
                           # data_universal_gene_summary = universal_gene_summary,
                           # var = "gene" #Variable to determine type of dependency table for cell lines
) {
  make_dep_table_raw <- function() {
    if(input$type == "gene") {
      data_universal_achilles_long <-
        get_data_object(object_names = input$content,
                        dataset_name = "universal_achilles_long",
                        pivotwider = TRUE) %>%
        dplyr::mutate(across(contains(c("score")), as.numeric)) %>%
        dplyr::select(-data_set)

      data_cell_expression_names <- get_content(object_name = "cell_expression_names", dataset = TRUE)

      table_data <-
        data_universal_achilles_long %>%
        tidyr::pivot_wider(names_from = id, values_from = dep_score) %>%
        dplyr::left_join(data_cell_expression_names, by = "depmap_id") %>%
        dplyr::select(cell_line, lineage, contains(input$content)) %>%
        dplyr::rename("Cell Line" = "cell_line", "Lineage" = "lineage")
      # } else if(input$type == "compound") {
      #   table_data <-
      #     data_universal_prism_long %>%
      #     dplyr::filter_all(dplyr::any_vars(name %in% input$content)) %>%
      #     tidyr::pivot_wider(names_from = name, values_from = log2fc) %>%
      #     dplyr::left_join(data_cell_expression_names, by = c("x1" = "X1")) %>%
      #     dplyr::select(cell_line, lineage, contains(input$content)) %>%
      #     dplyr::rename("Cell Line" = "cell_line", "Lineage" = "lineage")
      # } else {
      #   if(var == "gene"){
      #     table_data <-
      #       data_universal_achilles_long %>%
      #       dplyr::left_join(data_cell_expression_names, by = "X1") %>%
      #       dplyr::filter(cell_line %in% input$content) %>%
      #       dplyr::select(-X1, -lineage, -lineage_subtype) %>%
      #       tidyr::pivot_wider(names_from = cell_line, values_from = dep_score) %>%
      #       dplyr::left_join(data_universal_gene_summary, by = c("gene" = "approved_symbol")) %>%
      #       dplyr::select(gene, approved_name, contains(input$content)) %>%
      #       dplyr::mutate(unique_essential = gene %in% unique_essential_genes$gene, common_essential = gene %in% common_essentials$gene)
      #   }
      #   else if (var == "drug"){
      #     table_data <-
      #       data_universal_prism_long %>%
      #       dplyr::left_join(data_cell_expression_names, by = c("x1" = "X1")) %>%
      #       dplyr::left_join(data_compound_prism_names, by = "name") %>%
      #       dplyr::filter(cell_line %in% input$content) %>%
      #       tidyr::pivot_wider(names_from = cell_line, values_from = log2fc) %>%
      #       dplyr::rename("log2fc" = as.character(input$content)) %>%
      #       dplyr::select(name, moa, log2fc, dplyr::contains(input$content)) %>%
      #       dplyr::mutate(unique_toxic = name %in% prism_unique_toxic$name)
      # }
    }
    table_complete <-
      table_data %>%
      dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
      dplyr::arrange(.[[3]])

    return(table_complete)
  }
  #error handling
  tryCatch(make_dep_table_raw(),
           error = function(e){
             return(dplyr::tibble(`Cell Line` = NA,
                                  Lineage = NA,
                                  Gene = NA) %>%
                      dplyr::slice(0))
           })
}

## co-essentiality -----
#' Top Table
#'
#' \code{make_top_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a top Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a top Table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_top_table(input = list(type = 'gene', content = 'ROCK1'))
#' make_top_table(input = list(type = 'gene', content = c('ROCK1', 'ROCK2')))
#' \dontrun{
#' make_top_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_top_table <- function(input = list()) {
  make_top_table_raw <- function() {
    if(input$type == "gene") {
      table_data <-
        get_data_object(object_names = input$content,
                        dataset_name = "gene_master_top_table",
                        pivotwider = TRUE) %>%
        dplyr::mutate(dplyr::across(dplyr::contains(c("score", "r2", "concept")), as.numeric)) %>%
        dplyr::select(-data_set) %>%
        dplyr::arrange(dplyr::desc(r2))
      # } else if(input$type == "cell") {
      #   table_data <-
      #     data_master_top_table_cell %>%
      #     dplyr::filter_all(dplyr::any_vars(fav_cell %in% input$content))
    } else {
      stop("delcare your type")
    }
    return(table_data)
  }
  #error handling
  tryCatch(make_top_table_raw(),
           error = function(e){
             message(e)
           })
}

#' Bottom Table
#'
#' \code{make_bottom_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a bottom Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a bottom Table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_bottom_table(input = list(type = 'gene', content = 'ROCK1'))
#' make_bottom_table(input = list(type = 'gene', content = c('ROCK1', 'ROCK2')))
#' \dontrun{
#' make_bottom_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_bottom_table <- function(input = list(),
                              gls = FALSE,
                              ...) {
  make_bottom_table_raw <- function() {
    data_achilles_lower <- get_stats(data_set = "achilles", var = "lower")
    if(input$type == "gene") {
      table_data <-
        get_data_object(object_names = input$content,
                        dataset_name = "gene_master_bottom_table",
                        pivotwider = TRUE) %>%
        dplyr::mutate(dplyr::across(dplyr::contains(c("score", "r2", "concept")), as.numeric)) %>%
        dplyr::select(-data_set)
      # } else if(input$type == "cell") {
      #   table_data <-
      #     data_master_bottom_table_cell %>%
      #     dplyr::filter_all(dplyr::any_vars(fav_cell %in% input$content))
    } else {
      stop("delcare your type")
    }
    if(nrow(table_data) != 0) {
      if(gls){
        table_complete <-
          table_data %>%
          tidyr::unnest(data) %>%
          dplyr::filter(GLSpvalue < 0.05 & !duplicated(gene)) %>%
          dplyr::arrange(-dplyr::desc(GLSpvalue))
      } else {
        table_complete <-
          table_data %>%
          dplyr::filter(r2 < data_achilles_lower & !duplicated(gene)) %>%
          dplyr::arrange(dplyr::desc(r2))
      }
    } else {
      table_complete <-
        table_data
    }
    return(table_complete)
  }
  #error handling
  tryCatch(make_bottom_table_raw(),
           error = function(e){
             message(e)
           })
}

##censor-----
#' Censor
#'
#' Censor is used to remove genes from similarity table that are garbage (too many associations). It does this by filtering out genes above (greater_than) a provided threshold.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_censor_table(input = list(type = 'gene', content = 'ROCK1'))
#' make_censor_table(input = list(type = 'gene', content = 'ROCK1'), choice = TRUE)
#' make_censor_table(input = list(type = 'gene', content = 'ROCK1'), choice = TRUE, greater_than = 100)
#' \dontrun{
#' make_censor_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_censor_table <- function(input = list(),
                              censor = FALSE,
                              greater_than = 300) {
  table_data <-
    make_top_table(input = input)
  if(censor == TRUE){
    gene_censor_data <-
      get_data_object(object_names = input$content,
                      dataset_name = "gene_censor",
                      pivotwider = TRUE) %>%
      dplyr::mutate(across(contains(c("num_sim")), as.numeric)) %>%
      dplyr::select(-data_set)

    genes_to_censor <-
      gene_censor_data %>%
      dplyr::filter(num_sim > greater_than) %>% #get the genes that have too many associations
      dplyr::pull(id)

    censored_table <-
      table_data %>%
      dplyr::filter(!gene %in% genes_to_censor) #filter NOT in (out) the too high list

    return(censored_table)
  }
  return(table_data)
}

##enrichments-----
#' Enrichment Top Table
#'
#' \code{make_enrichment_top} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a enrichment top table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a enrichment top table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_gene_dependency_enrichment_table(input = list(type = 'gene', content = 'ROCK1'))
#' \dontrun{
#' make_gene_dependency_enrichment_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_gene_dependency_enrichment_table <- function(input = list()) {
  make_gene_dependency_enrichment_table_raw <- function() {
    gene_dependency_enrichment <-
      ddh::get_data_object(object_names = input$content,
                           dataset_name = "gene_dependency_enrichment",
                           pivotwider = TRUE) %>%
      dplyr::mutate(across(contains(c("NGenes", "PValue", "FDR")), as.numeric)) %>%
      dplyr::select(-data_set) %>%

      dplyr::mutate(`Gene Set` = gsub("_.*", "", GeneSet)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(GeneSet = gsub(paste0(`Gene Set`, "_"), "", GeneSet)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(GeneSet = stringr::str_replace_all(GeneSet, "_", " ")) %>%
      dplyr::mutate(GeneSet = stringr::str_to_sentence(GeneSet)) %>%
      dplyr::rename("Query" = "id", "Pathway" = "GeneSet", "# Genes" = "NGenes", "P-value" = "PValue")
    return(gene_dependency_enrichment)
  }
  #error handling
  tryCatch(make_gene_dependency_enrichment_table_raw(),
           error = function(e){
             message(e)
           })
}

#' Molecular Features Segments Table
#'
#' This is a table function that takes a gene name and returns a molecular features segments table
#'
#' @param input Expecting a list containing type and content variable.
#'
#' @return If no error, then returns a molecular features segments table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_molecular_features_segments_table(input = list(type = 'gene', content = 'ROCK1'))
#' make_molecular_features_segments_table(input = list(type = 'gene', content = c('ROCK1', 'ROCK2')))
#' \dontrun{
#' make_molecular_features_segments_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_molecular_features_segments_table <- function(input = list(),
                                                   ...) {
  make_molecular_features_segments_table_raw <- function() {
    gene_molecular_features_hits <-
      ddh::get_data_object(object_names = input$content,
                           dataset_name = "gene_molecular_features_segments",
                           pivotwider = TRUE) %>%
      dplyr::mutate(depscore = as.numeric(depscore)) %>%
      dplyr::select(-data_set) %>%
      dplyr::rename(Query = id)
    return(gene_molecular_features_hits)
  }
  #error handling
  tryCatch(make_molecular_features_segments_table_raw(),
           error = function(e){
             message(e)
           })
}

#' Molecular Features Table
#'
#' This is a table function that takes a gene name and returns a molecular features table
#'
#' @param input Expecting a list containing type and content variable.
#'
#' @return If no error, then returns a molecular features table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_molecular_features_table(input = list(type = 'gene', content = 'ROCK1'))
#' make_molecular_features_table(input = list(type = 'gene', content = c('ROCK1', 'ROCK2')))
#' \dontrun{
#' make_molecular_features_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_molecular_features_table <- function(input = list(),
                                          ...) {
  make_molecular_features_table_raw <- function() {
    gene_molecular_features_hits <-
      ddh::get_data_object(object_names = input$content,
                           dataset_name = "gene_molecular_features_top",
                           pivotwider = TRUE) %>%
      dplyr::mutate(dplyr::across(dplyr::contains(c("logFC", "pval", "adjPval")), as.numeric)) %>%
      dplyr::mutate_if(is.numeric, ~ signif(., digits = 3)) %>%
      dplyr::select(-data_set) %>%
      dplyr::rename(Query = id, Feature = feature, `P-value` = pval, FDR = adjPval)
    return(gene_molecular_features_hits)
  }
  #error handling
  tryCatch(make_molecular_features_table_raw(),
           error = function(e){
             message(e)
           })
}

#' Molecular Features Pathways Table
#'
#' This is a table function that takes a gene name and returns a molecular features pathways table
#'
#' @param input Expecting a list containing type and content variable.
#'
#' @return If no error, then returns a molecular features pathways table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_molecular_features_pathways_table(input = list(type = 'gene', content = 'ROCK1'))
#' make_molecular_features_pathways_table(input = list(type = 'gene', content = c('ROCK1', 'ROCK2')))
#' \dontrun{
#' make_molecular_features_pathways_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_molecular_features_pathways_table <- function(input = list(),
                                                   ...) {
  make_molecular_features_pathways_table_raw <- function() {
    gene_molecular_features_hits <-
      ddh::get_data_object(object_names = input$content,
                           dataset_name = "gene_molecular_features_pathways_top",
                           pivotwider = TRUE) %>%
      dplyr::mutate(dplyr::across(dplyr::contains(c("pval", "adjPval")), as.numeric)) %>%
      dplyr::mutate_if(is.numeric, ~ signif(., digits = 3)) %>%
      dplyr::select(-data_set) %>%
      dplyr::rename(Query = id, Pathway = GeneSet, `P-value` = pval, FDR = adjPval) %>%
      dplyr::mutate(`Gene Set` = gsub("_.*", "", Pathway)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Pathway = gsub(paste0(`Gene Set`, "_"), "", Pathway)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Pathway = stringr::str_replace_all(Pathway, "_", " ")) %>%
      dplyr::mutate(Pathway = stringr::str_to_sentence(Pathway))

    return(gene_molecular_features_hits)
  }
  #error handling
  tryCatch(make_molecular_features_pathways_table_raw(),
           error = function(e){
             message(e)
           })
}

#' Gene-Pathway CCA Table
#'
#' This is a table function that takes a gene name and returns a gene-pathway table
#'
#' @param input Expecting a list containing type and content variable.
#'
#' @return If no error, then returns a gene-pathway table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_cca_genes_table(input = list(type = 'gene', content = 'ROCK1'))
#' make_cca_genes_table(input = list(type = 'gene', content = c('ROCK1', 'ROCK2')))
#' \dontrun{
#' make_cca_genes_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_cca_genes_table <- function(input = list(),
                                 gene_set = NULL,
                                 ...) {
  make_cca_genes_table_raw <- function() {
    gene_pathway_hits <-
      ddh::get_data_object(object_names = input$content,
                           dataset_name = "gene_cca_pathway",
                           pivotwider = TRUE) %>%
      dplyr::mutate(dplyr::across(QuerySize:CCS, as.numeric)) %>%
      dplyr::select(-data_set) %>%
      dplyr::mutate(`Gene Set` = gsub("_.*", "", Pathway)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Pathway = gsub(paste0(`Gene Set`, "_"), "", Pathway)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Pathway = stringr::str_replace_all(Pathway, "_", " ")) %>%
      dplyr::mutate(Pathway = stringr::str_to_sentence(Pathway)) %>%
      dplyr::select(Query = id, Pathway, "# Genes" = PathwaySize,
                    "Explained Variance" = PathwayVarExp, CC = CC1, "Gene Set")

    if(!is.null(gene_set)) {
      if (!any(gene_set %in% gene_pathway_hits$`Gene Set`)) stop ("Gene set not found.")
      gene_pathway_hits <- gene_pathway_hits %>%
        dplyr::filter(`Gene Set` %in% gene_set)
    }

    return(gene_pathway_hits)
  }
  #error handling
  tryCatch(make_cca_genes_table_raw(),
           error = function(e){
             message(e)
           })
}

#' Pathway-Pathway CCA Table
#'
#' This is a table function that takes a pathway name and returns a pathway-pathway table
#'
#' @param input Expecting a list containing type and content variable.
#'
#' @return If no error, then returns a pathway-pathway table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_cca_pathway_table(input = list(type = 'pathway', content = 'GOBP_FIBROBLAST_GROWTH_FACTOR_PRODUCTION'))
#' \dontrun{
#' make_cca_pathway_table(input = list(type = 'pathway', content = 'GOBP_FIBROBLAST_GROWTH_FACTOR_PRODUCTION'))
#' }
make_cca_pathway_table <- function(input = list(),
                                   gene_set = NULL,
                                   ...) {
  make_cca_pathway_table_raw <- function() {
    pathway_pathway_hits <-
      ddh::get_data_object(object_names = input$content,
                           dataset_name = "pathway_cca_pathway",
                           pivotwider = TRUE) %>%
      dplyr::mutate(dplyr::across(QuerySize:CCS, as.numeric)) %>%
      dplyr::select(-data_set) %>%
      dplyr::mutate(`Gene Set` = gsub("_.*", "", Pathway)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Pathway = gsub(paste0(`Gene Set`, "_"), "", Pathway)) %>%
      dplyr::ungroup() %>%
      dplyr::rename(Query = id) %>%
      dplyr::mutate(Query = stringr::str_replace_all(Query, "_", " "),
                    Pathway = stringr::str_replace_all(Pathway, "_", " ")) %>%
      dplyr::mutate(Query = stringr::str_to_sentence(Query),
                    Pathway = stringr::str_to_sentence(Pathway))

    if(!is.null(gene_set)) {
      if (!any(gene_set %in% pathway_pathway_hits$`Gene Set`)) stop ("Gene set not found.")
      pathway_pathway_hits <- pathway_pathway_hits %>%
        dplyr::filter(`Gene Set` %in% gene_set)
    }

    return(pathway_pathway_hits)
  }
  #error handling
  tryCatch(make_cca_pathway_table_raw(),
           error = function(e){
             message(e)
           })
}

# CELL SUMMARY TABLE ------
#' Cell Summary Table
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_cell_line_table(input = list(content = 'HEPG2'))
make_cell_line_table <- function(#data_cell_expression_meta = cell_expression_meta,
  input = list()) {

  make_cell_summary_raw <- function() {

    # queried_lineage <-
    #   data_cell_expression_meta %>%
    #   dplyr::filter(cell_line %in% input$content) %>%
    #   dplyr::pull(lineage) %>%
    #   unique()
    #
    # cell_summary_table <-
    #   data_cell_expression_meta %>%
    #   dplyr::filter(!cell_line %in% input$content) %>%
    #   dplyr::filter(lineage %in% queried_lineage) %>%
    #   dplyr::select(`Cell Line` = cell_line,
    #                 Lineage = lineage,
    #                 `Lineage Subtype` = lineage_subtype,
    #                 Age = age,
    #                 Sex = sex)
    #
    # return(cell_summary_table)
  }

  # error handling
  tryCatch(make_cell_summary_raw(),
           error = function(e){
             message(e)
           })

}

# CELL SIMILARITY TABLE -----
#' Cell Similarity Table
#'
#' \code{make_cell_sim_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a cell sim Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a cell sim Table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_cell_sim_table(input = list(type = "cell", content = "HEL"))
#' make_cell_sim_table(input = list(type = "cell", content = c("HEL", "LS513")))
#' \dontrun{
#' make_cell_sim_table(input = list(type = "cell", content = "HEL"))
#' }
make_cell_sim_table <- function(input = list(),
                                #data_cell_dependency_sim = cell_dependency_sim,
                                #data_cell_expression_sim = cell_expression_sim,
                                similarity = "dependency",
                                bonferroni_cutoff = 0.05) {
  make_cell_sim_table_raw <- function() {
    glue::glue('{similarity} and {bonferroni_cutoff}')
    # if(similarity == "dependency") {
    #   cell_sims <- data_cell_dependency_sim
    # } else if(similarity == "expression") {
    #   cell_sims <- data_cell_expression_sim
    # }
    #
    # cell_sim_table <-
    #   cell_sims %>%
    #   dplyr::filter(cell1_name %in% input$content | cell2_name %in% input$content) %>%
    #   dplyr::mutate(sex = ifelse(cell1_name %in% input$content, sex.y, sex.x),
    #                 age = ifelse(cell1_name %in% input$content, age.y, age.x),
    #                 lineage = ifelse(cell1_name %in% input$content, lineage.y, lineage.x),
    #                 lineage_subtype = ifelse(cell1_name %in% input$content, lineage_subtype.y, lineage_subtype.x),
    #                 status = ifelse(cell1_name %in% input$content, primary_or_metastasis.y, primary_or_metastasis.x)
    #   ) %>%
    #   dplyr::select(cell1_name, cell2_name, lineage, coef, pval, fdr, bonferroni, lineage_subtype, sex, age, status)
    #
    # # Filter significant
    # cell_sim_table <-
    #   cell_sim_table %>%
    #   dplyr::filter(bonferroni < bonferroni_cutoff)
    #
    # # Swap cols (based on query)
    # for(i in 1:nrow(cell_sim_table)) {
    #   if(cell_sim_table$cell2_name[i] %in% input$content & !(cell_sim_table$cell1_name[i] %in% input$content)) {
    #     cell1 <- cell_sim_table$cell1_name[i]
    #     cell2 <- cell_sim_table$cell2_name[i]
    #
    #     cell_sim_table$cell2_name[i] <- cell1
    #     cell_sim_table$cell1_name[i] <- cell2
    #   }
    # }
    #
    # # top table
    # top_table <-
    #   cell_sim_table %>%
    #   dplyr::filter(coef > 0)
    #
    # # bottom table
    # bottom_table <-
    #   cell_sim_table %>%
    #   dplyr::filter(coef < 0)
    #
    # return(list(top_table = dplyr::as_tibble(top_table),
    #             bottom_table = dplyr::as_tibble(bottom_table)))
  }

  # error handling
  tryCatch(make_cell_sim_table_raw(),
           error = function(e){
             message(e)
           })
}

# DRUG TABLES -----
#' Drug Genes Cor Table
#'
#' \code{make_drug_genes_cor_table} takes a drug name and returns a table of genes
#'
#' This is a table function that takes a drug name and returns a table of genes
#'
#' @importFrom magrittr %>%
#'
#' @export
make_drug_genes_cor_table <- function(input = list()
                                      #data_compound_genes_cor_table = compound_genes_cor_table,
                                      #drug
) {
  make_drug_genes_cor_table_raw <- function() {
    # unnested_table <-
    #   data_compound_genes_cor_table %>%
    #   dplyr::filter_all(dplyr::any_vars(fav_drug %in% input$content)) %>%
    #   tidyr::unnest(data) %>%
    #   dplyr::ungroup() %>%
    #   dplyr::arrange(dplyr::desc(r2))
    # return(unnested_table)
  }
  tryCatch(make_drug_genes_cor_table_raw(),
           error = function(e){
             message(e)
           })
}

#' Gene Drugs Cor Table
#'
#' This is a table function that takes a gene name and returns a gene drugs cor Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a gene drugs cor Table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_gene_drugs_cor_table(input = list(type = 'gene', content = 'ACOT4'))
#' \dontrun{
#' make_gene_drugs_cor_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_gene_drugs_cor_table <- function(input = list(),
                                      ...) {
  make_gene_drugs_cor_table_raw <- function() {
    gene_drugs_cor_table <-
      get_data_object(object_names = input$content,
                      dataset_name = "gene_drugs_cor_table",
                      pivotwider = TRUE)

    if(nrow(gene_drugs_cor_table) > 0) {
      gene_drugs_cor_table <-
        gene_drugs_cor_table%>%
        dplyr::mutate(across(contains(c("z_score", "r2")), as.numeric)) %>%
        dplyr::select(-data_set) %>%
        dplyr::arrange(dplyr::desc(r2)) %>%
        dplyr::rename(Query = id, Drug = drug, MOA = moa, `Z-score` = z_score, R2 = r2)
    }
    return(gene_drugs_cor_table)
  }
  tryCatch(make_gene_drugs_cor_table_raw(),
           error = function(e){
             message(e)
           })
}

#' Drug Genes Table
#'
#' \code{make_drug_genes_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a drug genes Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a drug genes Table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_drug_genes_table(drug = "aspirin")
#' make_drug_genes_table(drug = "ibuprofen")
#' make_drug_genes_table(drug = c("aspirin", "ibuprofen"))
#' \dontrun{
#' make_drug_genes_table(drug = "aspirin")
#' }
make_drug_genes_table <- function(input = list()
                                  #data_compound_genes_table = compound_genes_table,
                                  #drug
) {
  make_drug_genes_table_raw <- function() {
    # unnested_table <-
    #   data_compound_genes_table %>%
    #   dplyr::filter_all(dplyr::any_vars(fav_drug %in% input$content)) %>%
    #   tidyr::unnest(data) %>%
    #   dplyr::ungroup() %>%
    #   dplyr::arrange(fav_gene)
    # return(unnested_table)
    NULL
  }
  #error handling
  tryCatch(make_drug_genes_table_raw(),
           error = function(e){
             message(e)
           })
}

#' Gene Drugs Table
#'
#' This is a table function that takes a gene name and returns a drug genes Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a drug genes Table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_gene_drugs_table(input = list(type = "gene", content = "ROCK1"))
#' \dontrun{
#' make_gene_drugs_table(input = list(type = "gene", content = "ROCK1"))
#' }
make_gene_drugs_table <- function(input = list()) {
  make_gene_drugs_table_raw <- function() {
    gene_drugs_table <-
      ddh::get_data_object(object_names = input$content,
                           dataset_name = "gene_drugs_table",
                           pivotwider = TRUE) %>%
      dplyr::select(all_of(c("id", "fav_drug", "moa"))) %>%
      dplyr::arrange(fav_drug)
    return(gene_drugs_table)
  }
  #error handling
  tryCatch(make_gene_drugs_table_raw(),
           error = function(e){
             message(e)
           })
}

#' Cell Drugs Table
#'
#' @importFrom magrittr %>%
#'
#' @export
make_cell_drugs_table <- function(#data_universal_prism_long = universal_prism_long,
  #data_cell_expression_meta = cell_expression_meta,
  input = list()) {
  make_cell_drugs_table_raw <- function() {
    # cell_table <-
    #   data_universal_prism_long %>%
    #   dplyr::rename(X1 = 1) %>%
    #   dplyr::left_join(data_cell_expression_meta, by = "X1") %>%
    #   dplyr::select(cell_line, name, log2fc) %>%
    #   dplyr::filter(cell_line %in% input$content) %>%
    #   dplyr::arrange(dplyr::desc(abs(log2fc)))
    # return(cell_table)
  }
  #error handling
  tryCatch(make_cell_drugs_table_raw(),
           error = function(e){
             message(e)
           })
}

# METABOLITE TABLES -----
#' Metabolite Table
#'
#' \code{make_metabolite_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a metabolite Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a metabolite Table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_metabolite_table(input = list(type = 'gene', content = 'ROCK1'))
#' make_metabolite_table(input = list(type = 'gene', content = 'ROCK1'), collapse = TRUE)
#' make_metabolite_table(input = list(type = "gene", content = "ROCK2"))
#' make_metabolite_table(input = list(type = "gene", content = c("ROCK1", "ROCK2")))
#' make_metabolite_table(input = list(type = "compound", content = "NAD"))
#' make_metabolite_table(input = list(type = "cell", content = "HEPG2"))
#' \dontrun{
#' make_metabolite_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_metabolite_table <- function(input = list(),
                                  collapse = FALSE) {
  make_metabolite_table_raw <- function() {
    var <- dplyr::if_else(collapse == TRUE, "metabolite_collapsed", "metabolite")
    if(input$type == "gene") {
      table_complete <-
        get_data_object(object_names = input$content,
                        dataset_name = "compound_hmdb_proteins",
                        pivotwider = TRUE) %>%
        dplyr::select(all_of(c("id", var))) %>%
        tidyr::drop_na(all_of(var))
    } else if(input$type == "compound") {
      # table_complete <-
      #   data_compound_hmdb_metabolites %>% #metabolite centric df
      #   dplyr::filter_all(dplyr::any_vars(fav_metabolite %in% input$content)) %>%
      #   tidyr::unnest(data) %>%
      #   dplyr::ungroup() %>%
      #   dplyr::select(metabolite_name, gene_name, gene_accession, metabolite_accession) %>%
      #   dplyr::arrange(metabolite_name, gene_name)
    } else if(input$type == "cell") {
      # table_complete <-
      #   data_cell_metabolites %>%
      #   dplyr::left_join(data_cell_expression_names, by = c("DepMap_ID" = "X1")) %>%
      #   dplyr::filter(cell_line %in% input$content) %>%
      #   dplyr::select(-CCLE_ID, -DepMap_ID, -lineage, -lineage_subtype) %>%
      #   tidyr::pivot_longer(cols = !cell_line, names_to = "metabolite") %>%
      #   dplyr::arrange(-value)
    } else {
      stop("declare your type") }
    return(table_complete)
  }
  #error handling
  tryCatch(make_metabolite_table_raw(),
           error = function(e){
             message(e)
           })
}
