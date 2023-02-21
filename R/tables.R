
# PATHWAY TABLES -----
#' Pathway List Table
#'
#' \code{make_pathway_list} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a pathway list table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a pathway list table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_pathway_list(input = list(type = 'gene', content = 'ROCK1'))
#' \dontrun{
#' make_pathway_list(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_pathway_list <- function(data_gene_pathways = gene_pathways,
                              input = list()) { #makes a subtable of gene_pathways that contains your gene query
  make_pathway_list_raw <- function() {
    present <- function(list, query){ #is the gene present?
      y <- unlist(list, use.names = FALSE)
      any(stringr::str_detect(y, query))
    }
    filtered_table <-
      data_gene_pathways %>%
      dplyr::filter(purrr::map_lgl(data_gene_pathways$data, ~present(list = .x, query = input$content)) == TRUE)
    return(filtered_table)
  }
  #error handling
  tryCatch(make_pathway_list_raw(),
           error = function(e){
             #make empty table equivalent to returned table
             return(data_gene_pathways %>%
                      dplyr::slice(0)
             )
           })
}

#' Pathway Genes Table
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_pathway_genes(go_id = "1902965")
make_pathway_genes <- function(data_gene_pathways = gene_pathways,
                               data_universal_gene_summary = universal_gene_summary,
                               go_id) {
  make_pathway_genes_raw <- function() {
    pathway_table <-
      data_gene_pathways %>%
      dplyr::filter(go %in% go_id) %>%
      tidyr::unnest(data) %>%
      dplyr::left_join(data_universal_gene_summary, by = c("gene" = "approved_symbol")) %>%
      dplyr::select(gene, approved_name, aka)
    return(pathway_table)
  }
  #error handling
  tryCatch(make_pathway_genes_raw(),
           error = function(e){
             #make empty table equivalent to returned table
             return(data_gene_pathways %>%
                      tidyr::unnest(data) %>%
                      dplyr::left_join(data_universal_gene_summary, by = c("gene" = "approved_symbol")) %>%
                      dplyr::select(gene, approved_name, aka) %>%
                      dplyr::slice(0))
           })
}

#' Compound Table
#'
#' \code{make_compound_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a compound Table
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
make_compound_table <- function(data_compound_prism_cor_nest = compound_prism_cor_nest,
                                data_compound_prism_names = compound_prism_names,
                                input = list(),
                                data_universal_stats_summary = universal_stats_summary,
                                top = TRUE) {
  make_compound_table_raw <- function() {
    prism_cor_upper <- get_stats(data_set = "prism", var = "upper")
    prism_cor_lower <- get_stats(data_set = "prism", var = "lower")

    table_complete <-
      data_compound_prism_cor_nest %>%
      dplyr::filter_all(dplyr::any_vars(fav_drug %in% input$content)) %>%
      tidyr::unnest("data") %>%
      dplyr::ungroup() %>%
      {if (top == TRUE) dplyr::filter(., r2 > prism_cor_upper) else dplyr::filter(., r2 < prism_cor_lower)} %>%
      dplyr::arrange(dplyr::desc(r2)) %>%
      dplyr::left_join(data_compound_prism_names, by = "name") %>%
      dplyr::select(1:4)
    return(table_complete)
  }
  #error handling
  tryCatch(make_compound_table_raw(),
           error = function(e){
             #make empty table equivalent to returned table
             return(data_prism_cor_nest %>%
                      tidyr::unnest("data") %>%
                      dplyr::ungroup() %>%
                      dplyr::left_join(data_compound_prism_names, by = "name") %>%
                      dplyr::select(1:4) %>%
                      dplyr::slice(0))
           }
  )
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
      dplyr::rename(Gene = id,
                    Reliability = reliability,
                    Location = main_location,
                    Expression = value) %>%
      dplyr::arrange(Expression)
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
      dplyr::select(-name)

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

## PROTEIN CLUSTER TABLE -----------------------------------------------
#' Clustering Table
#'
#' \code{make_clustering_table} returns an image of ...
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
#' make_clustering_table(input = list(type = 'gene', content = 'ROCK1'))
#' \dontrun{
#' make_clustering_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_clustering_table <- function(input = list(),
                                  cluster = FALSE,
                                  show_signature = FALSE) {
  make_clustering_table_raw <- function() {

    data_gene_signature_clusters <-
      get_data_object(object_names = input$content,
                    dataset_name = "gene_signature_clusters",
                    pivotwider = TRUE) %>%
      dplyr::mutate(across(contains(c("X", "clust", "member_prob")), as.numeric))

    data_gene_signatures <-
      get_data_object(object_names = input$content,
                      dataset_name = "gene_signatures",
                      pivotwider = TRUE) %>%
      dplyr::mutate(across(A:Y, as.numeric))

    gene_signature_clusters <- get_content("gene_signature_clusters", dataset = TRUE) #allows all genes, coordinates, and cluster
    gene_signatures <- get_content("gene_signatures", dataset = TRUE) #all gene signatures

    #vec of clusters in query
    query_clust <-
      data_gene_signature_clusters %>%
      dplyr::pull(clust) %>%
      unique()

    if(show_signature == TRUE){cluster <- TRUE} #make sure you preserve clusters
    #table of all proteins in cluster containing: uniprot_id, gene_name, protein_name, clust, cluster_name
    cluster_table <-
      gene_signature_clusters %>%
      dplyr::filter(clust %in% query_clust) %>%
      dplyr::select(gene_name = id, protein_name, description) %>%
      #if cluster is true, then all proteins; if false, then just the query
      {if(cluster == FALSE) dplyr::filter(., gene_name %in% input$content) else .}

    # if show_sig is true, then calculate the signature of the cluster on the fly
    if(show_signature) {
      cluster_table <-
        cluster_table %>%
        dplyr::left_join(gene_signatures %>%
                           dplyr::select(gene_name, A:U) %>%
                           dplyr::mutate_at(dplyr::vars(A:U), ~ round(., 2)),
                         by = "gene_name")
    }

    if(nrow(cluster_table) == 0) {
      stop("Unable to cluster this protein by its amino acid sequence so far...")
    }
    return(dplyr::as_tibble(cluster_table))
  }

  #error handling
  tryCatch(make_clustering_table_raw(),
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

    gene_signature_clusters <- get_content("gene_signature_clusters", dataset = TRUE) #allows all genes, coordinates, and cluster
    gene_signature_cluster_enrichment <- get_content("gene_signature_cluster_enrichment", dataset = TRUE)

    if(filter_noise) {
      gene_signature_clusters <-
        gene_signature_clusters %>%
        dplyr::filter(clust != 0)
    }

    query_clust <-
      gene_signature_clusters %>%
      dplyr::filter(id %in% input$content) %>%
      dplyr::pull(clust) %>%
      unique()

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
make_structure3d_table <- function(data_gene_uniprot_pdb_table = gene_uniprot_pdb_table,
                                   data_universal_proteins = universal_proteins,
                                   input = list()) {
  make_structure3d_table_raw <- function() {
    table_complete <-
      data_universal_proteins %>%
      dplyr::filter(gene_name %in% input$content) %>%
      dplyr::left_join(data_gene_uniprot_pdb_table, by = c("uniprot_id" = "uniprot")) %>%
      tidyr::unnest(data) %>%
      dplyr::select(gene_name, uniprot_id, pdb, title, doi, organism, expression_system)

    return(table_complete)
  }

  #error handling
  tryCatch(make_structure3d_table_raw(),
           error = function(e){
             return(data_compound_hmdb_proteins %>%
                      dplyr::left_join(data_gene_uniprot_pdb_table, by = c("uniprot_id" = "uniprot")) %>%
                      tidyr::unnest(data) %>%
                      dplyr::select(gene_name, uniprot_id, pdb, title, doi, organism, expression_system) %>%
                      dplyr::slice(0))
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
#' make_dep_table(input = list(type = 'compound', content = 'aspirin'))
#' make_dep_table(input = list(type = 'cell', content = 'HEPG2'))
#' \dontrun{
#' make_dep_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_dep_table <- function(data_universal_achilles_long = universal_achilles_long,
                           data_universal_prism_long = universal_prism_long,
                           data_compound_prism_names = compound_prism_names,
                           data_cell_expression_names = cell_expression_names,
                           data_universal_gene_summary = universal_gene_summary,
                           input = list(),
                           var = "gene") { #Variable to determine type of dependency table for cell lines
  make_dep_table_raw <- function() {
    if(input$type == "gene") {
      table_data <-
        data_universal_achilles_long %>%
        dplyr::filter_all(dplyr::any_vars(gene %in% input$content)) %>%
        tidyr::pivot_wider(names_from = gene, values_from = dep_score) %>%
        dplyr::left_join(data_cell_expression_names, by = "X1") %>%
        dplyr::select(cell_line, lineage, contains(input$content)) %>%
        dplyr::rename("Cell Line" = "cell_line", "Lineage" = "lineage")
    } else if(input$type == "compound") {
      table_data <-
        data_universal_prism_long %>%
        dplyr::filter_all(dplyr::any_vars(name %in% input$content)) %>%
        tidyr::pivot_wider(names_from = name, values_from = log2fc) %>%
        dplyr::left_join(data_cell_expression_names, by = c("x1" = "X1")) %>%
        dplyr::select(cell_line, lineage, contains(input$content)) %>%
        dplyr::rename("Cell Line" = "cell_line", "Lineage" = "lineage")
    } else {
      if(var == "gene"){
        table_data <-
          data_universal_achilles_long %>%
          dplyr::left_join(data_cell_expression_names, by = "X1") %>%
          dplyr::filter(cell_line %in% input$content) %>%
          dplyr::select(-X1, -lineage, -lineage_subtype) %>%
          tidyr::pivot_wider(names_from = cell_line, values_from = dep_score) %>%
          dplyr::left_join(data_universal_gene_summary, by = c("gene" = "approved_symbol")) %>%
          dplyr::select(gene, approved_name, contains(input$content)) %>%
          dplyr::mutate(unique_essential = gene %in% unique_essential_genes$gene, common_essential = gene %in% common_essentials$gene)
      }
      else if (var == "drug"){
        table_data <-
          data_universal_prism_long %>%
          dplyr::left_join(data_cell_expression_names, by = c("x1" = "X1")) %>%
          dplyr::left_join(data_compound_prism_names, by = "name") %>%
          dplyr::filter(cell_line %in% input$content) %>%
          tidyr::pivot_wider(names_from = cell_line, values_from = log2fc) %>%
          dplyr::rename("log2fc" = as.character(input$content)) %>%
          dplyr::select(name, moa, log2fc, dplyr::contains(input$content)) %>%
          dplyr::mutate(unique_toxic = name %in% prism_unique_toxic$name)
      }
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
#' \dontrun{
#' make_top_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_top_table <- function(data_gene_master_top_table = gene_master_top_table,
                           #data_master_top_table_cell = master_top_table_cell_line,
                           data_universal_stats_summary = universal_stats_summary,
                           gls = FALSE,
                           input = list()) {
  make_top_table_raw <- function() {
    data_achilles_upper <- get_stats(data_set = "achilles", var = "upper")
    if(input$type == "gene") {
      table_data <-
        data_gene_master_top_table %>%
        dplyr::filter_all(dplyr::any_vars(fav_gene %in% input$content))
      # } else if(input$type == "cell") {
      #   table_data <-
      #     data_master_top_table_cell %>%
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
          tidyr::unnest(data) %>%
          dplyr::filter(r2 > data_achilles_upper & !duplicated(gene)) %>%
          dplyr::arrange(dplyr::desc(r2))
      }
    } else {
      table_complete <-
        table_data
    }
    return(table_complete)
  }
  #error handling
  tryCatch(make_top_table_raw(),
           error = function(e){
             return(data_gene_master_top_table %>%
                      tidyr::unnest(data) %>%
                      dplyr::slice(0)
             )
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
#' \dontrun{
#' make_bottom_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_bottom_table <- function(data_gene_master_bottom_table = gene_master_bottom_table,
                              # data_master_bottom_table_cell = master_bottom_table_cell_line,
                              data_universal_stats_summary = universal_stats_summary,
                              gls = FALSE,
                              input = list()) {
  make_bottom_table_raw <- function() {
    data_achilles_lower <- get_stats(data_set = "achilles", var = "lower")
    if(input$type == "gene") {
      table_data <-
        data_gene_master_bottom_table %>%
        dplyr::filter_all(dplyr::any_vars(fav_gene %in% input$content))
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
          tidyr::unnest(data) %>%
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
             return(data_gene_master_bottom_table %>%
                      tidyr::unnest(data) %>%
                      dplyr::slice(0))
           })
}

#' Gene-Pathway and Pathway-Pathway Co-essentiality Table
#'
#' This is a table function that takes a gene/s or a pathway/s and returns its co-essential pathways and genes
#'
#' @param input Expecting a list containing content variable.
#' @param cutoff Absolute Pearson's correlation value to filter
#'
#' @return If no error, then returns a table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_gene_pathways_components(input = list(type = 'gene', content = 'ROCK1'))
#' \dontrun{
#' make_gene_pathways_components(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_gene_pathways_components <- function(data_gene_pathways_components = gene_pathways_components,
                                          input = list(),
                                          cutoff = NULL) {
  make_gene_pathways_components_raw <- function() {

    table_complete <- data_gene_pathways_components %>%
      dplyr::filter(feature1 %in% input$content | feature2 %in% input$content)

    if(is.null(cutoff)) {
      cutoff <- mean(abs(table_complete$pearson_corr)) # + sd(abs(table_complete$pearson_corr))
    }

    table_complete <- table_complete %>%
      dplyr::filter(abs(pearson_corr) > cutoff) %>%
      dplyr::mutate(swapped = FALSE) %>%
      dplyr::as_tibble()

    # Swap cols (based on query)
    for(i in 1:nrow(table_complete)) {
      if(table_complete$feature2[i] %in% input$content &
         !(table_complete$feature1[i] %in% input$content)) {
        ft1 <- table_complete$feature1[i]
        ft2 <- table_complete$feature2[i]

        table_complete$feature2[i] <- ft1
        table_complete$feature1[i] <- ft2
        table_complete$swapped[i] <- TRUE
      }
    }

    return(table_complete)
  }

  #error handling
  tryCatch(make_gene_pathways_components_raw(),
           error = function(e){
             return(data_gene_pathways_components %>%
                      dplyr::slice(0)
             )
           })
}

##censor-----
#' Censor
#'
#' Censor is used to remove genes from similarity table that are garbage (too many associations)
#'
#' @importFrom magrittr %>%
#'
#' @export
censor <- function(data_gene_master_top_table = gene_master_top_table,
                   data_gene_censor = gene_censor,
                   choice = FALSE,
                   greater_than) {
  if(choice == TRUE){
    data_gene_censor <-
      data_gene_censor %>%
      dplyr::filter(num_sim > greater_than) #get the genes that have too many associations

    censored_table <-
      data_gene_master_top_table %>%
      dplyr::filter(!gene %in% data_gene_censor$genes) #filter NOT in (out) the too high list

    return(censored_table)
  }
  return(data_gene_master_top_table)
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
#' make_enrichment_top(input = list(type = 'gene', content = 'ROCK1'))
#' \dontrun{
#' make_enrichment_top(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_enrichment_top <- function(data_gene_master_positive = gene_master_positive,
                                input = list()) {
  make_enrichment_top_raw <- function() {
    data_gene_master_positive %>%
      dplyr::filter_all(dplyr::any_vars(fav_gene %in% input$content)) %>%
      tidyr::unnest(data) %>%
      dplyr::select(fav_gene, enrichr, Term, Overlap, Adjusted.P.value, Combined.Score, Genes) %>%
      dplyr::arrange(Adjusted.P.value) %>%
      dplyr::rename("Query" = "fav_gene", "Gene Set" = "enrichr", "Gene List" = "Term", "Adjusted p-value" = "Adjusted.P.value", "Combined Score" = "Combined.Score") #"Overlap", "Genes"
  }
  #error handling
  tryCatch(make_enrichment_top_raw(),
           error = function(e){
             return(data_gene_master_positive %>%
                      tidyr::unnest(data) %>%
                      dplyr::select(fav_gene, enrichr, Term, Overlap, Adjusted.P.value, Combined.Score, Genes) %>%
                      dplyr::rename("Query" = "fav_gene",
                                    "Gene Set" = "enrichr",
                                    "Gene List" = "Term",
                                    "Adjusted p-value" = "Adjusted.P.value",
                                    "Combined Score" = "Combined.Score") %>%
                      dplyr::slice(0))
           })
}

#' Enrichment Bottom Table
#'
#' \code{make_enrichment_bottom} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a enrichment bottom table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a enrichment bottom table. If an error is thrown, then will return an empty table.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_enrichment_bottom(input = list(type = 'gene', content = 'ROCK1'))
#' \dontrun{
#' make_enrichment_bottom(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_enrichment_bottom <- function(data_gene_master_negative = gene_master_negative,
                                   input = list()) {
  make_enrichment_bottom_raw <- function() {
    data_gene_master_negative %>%
      dplyr::filter_all(dplyr::any_vars(fav_gene %in% input$content)) %>%
      tidyr::unnest(data) %>%
      dplyr::select(fav_gene, enrichr, Term, Overlap, Adjusted.P.value, Combined.Score, Genes) %>%
      dplyr::arrange(Adjusted.P.value) %>%
      dplyr::rename("Query" = "fav_gene", "Gene Set" = "enrichr", "Gene List" = "Term", "Adjusted p-value" = "Adjusted.P.value", "Combined Score" = "Combined.Score") #"Overlap", "Genes"
  }
  #error handling
  tryCatch(make_enrichment_bottom_raw(),
           error = function(e){
             return(data_gene_master_negative %>%
                      tidyr::unnest(data) %>%
                      dplyr::select(fav_gene, enrichr, Term, Overlap, Adjusted.P.value, Combined.Score, Genes) %>%
                      dplyr::rename("Query" = "fav_gene",
                                    "Gene Set" = "enrichr",
                                    "Gene List" = "Term",
                                    "Adjusted p-value" = "Adjusted.P.value",
                                    "Combined Score" = "Combined.Score") %>%
                      dplyr::slice(0))
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
make_cell_line_table <- function(data_cell_expression_meta = cell_expression_meta,
                                 input = list()) {

  make_cell_summary_raw <- function() {

    queried_lineage <-
      data_cell_expression_meta %>%
      dplyr::filter(cell_line %in% input$content) %>%
      dplyr::pull(lineage) %>%
      unique()

    cell_summary_table <-
      data_cell_expression_meta %>%
      dplyr::filter(!cell_line %in% input$content) %>%
      dplyr::filter(lineage %in% queried_lineage) %>%
      dplyr::select(`Cell Line` = cell_line,
                    Lineage = lineage,
                    `Lineage Subtype` = lineage_subtype,
                    Age = age,
                    Sex = sex)

    return(cell_summary_table)
  }

  # error handling
  tryCatch(make_cell_summary_raw(),
           error = function(e){
             return(dplyr::tibble(`Cell Line` = NA,
                                  Lineage = NA,
                                  `Lineage Subtype` = NA,
                                  Age = NA,
                                  Sex = NA) %>%
                      dplyr::slice(0))
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
make_cell_sim_table <- function(data_cell_dependency_sim = cell_dependency_sim,
                                data_cell_expression_sim = cell_expression_sim,
                                similarity = "dependency",
                                bonferroni_cutoff = 0.05,
                                input = list()) {
  make_cell_sim_table_raw <- function() {
    if(similarity == "dependency") {
      cell_sims <- data_cell_dependency_sim
    } else if(similarity == "expression") {
      cell_sims <- data_cell_expression_sim
    }

    cell_sim_table <-
      cell_sims %>%
      dplyr::filter(cell1_name %in% input$content | cell2_name %in% input$content) %>%
      dplyr::mutate(sex = ifelse(cell1_name %in% input$content, sex.y, sex.x),
                    age = ifelse(cell1_name %in% input$content, age.y, age.x),
                    lineage = ifelse(cell1_name %in% input$content, lineage.y, lineage.x),
                    lineage_subtype = ifelse(cell1_name %in% input$content, lineage_subtype.y, lineage_subtype.x),
                    status = ifelse(cell1_name %in% input$content, primary_or_metastasis.y, primary_or_metastasis.x)
      ) %>%
      dplyr::select(cell1_name, cell2_name, lineage, coef, pval, fdr, bonferroni, lineage_subtype, sex, age, status)

    # Filter significant
    cell_sim_table <-
      cell_sim_table %>%
      dplyr::filter(bonferroni < bonferroni_cutoff)

    # Swap cols (based on query)
    for(i in 1:nrow(cell_sim_table)) {
      if(cell_sim_table$cell2_name[i] %in% input$content & !(cell_sim_table$cell1_name[i] %in% input$content)) {
        cell1 <- cell_sim_table$cell1_name[i]
        cell2 <- cell_sim_table$cell2_name[i]

        cell_sim_table$cell2_name[i] <- cell1
        cell_sim_table$cell1_name[i] <- cell2
      }
    }

    # top table
    top_table <-
      cell_sim_table %>%
      dplyr::filter(coef > 0)

    # bottom table
    bottom_table <-
      cell_sim_table %>%
      dplyr::filter(coef < 0)

    return(list(top_table = dplyr::as_tibble(top_table),
                bottom_table = dplyr::as_tibble(bottom_table)))
  }

  # error handling
  tryCatch(make_cell_sim_table_raw(),
           error = function(e){
             return(dplyr::tibble(cell1_name = NA, cell2_name = NA,
                                  lineage = NA, coef = NA, pval = NA,
                                  fdr = NA, bonferroni = NA, lineage_subtype = NA,
                                  sex = NA, age = NA, status = NA) %>%
                      dplyr::slice(0))
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
make_drug_genes_cor_table <- function(data_compound_genes_cor_table = compound_genes_cor_table,
                                      drug) {
  make_drug_genes_cor_table_raw <- function() {
    unnested_table <-
      data_compound_genes_cor_table %>%
      dplyr::filter_all(dplyr::any_vars(fav_drug %in% drug)) %>%
      tidyr::unnest(data) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(r2))
    return(unnested_table)
  }
  tryCatch(make_drug_genes_cor_table_raw(),
           error = function(e){
             #make empty table equivalent to returned table
             return(data_compound_genes_cor_table %>%
                      tidyr::unnest(cols = c(data)) %>%
                      dplyr::slice(0))
           })
}

#' Gene Drugs Cor Table
#'
#' \code{make_gene_drugs_cor_table} returns an image of ...
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
#' make_gene_drugs_cor_table(input = list(type = 'gene', content = 'ROCK1'))
#' \dontrun{
#' make_gene_drugs_cor_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_gene_drugs_cor_table <- function(data_gene_drugs_cor_table = gene_drugs_cor_table,
                                      input = list()) {
  make_gene_drugs_cor_table_raw <- function() {
    unnested_table <-
      data_gene_drugs_cor_table %>%
      dplyr::filter_all(dplyr::any_vars(fav_gene %in% input$content)) %>%
      tidyr::unnest(data) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(r2))
    return(unnested_table)
  }
  tryCatch(make_gene_drugs_cor_table_raw(),
           error = function(e){
             #make empty table equivalent to returned table
             return(data_gene_drugs_cor_table %>%
                      tidyr::unnest(cols = c(data)) %>%
                      dplyr::slice(0))
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
make_drug_genes_table <- function(data_compound_genes_table = compound_genes_table,
                                  drug) {
  make_drug_genes_table_raw <- function() {
    unnested_table <-
      data_compound_genes_table %>%
      dplyr::filter_all(dplyr::any_vars(fav_drug %in% drug)) %>%
      tidyr::unnest(data) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(fav_gene)
    return(unnested_table)
  }
  #error handling
  tryCatch(make_drug_genes_table_raw(),
           error = function(e){
             #make empty table equivalent to returned table
             return(data_compound_genes_table %>%
                      tidyr::unnest(cols = c(data)) %>%
                      dplyr::slice(0))
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
make_gene_drugs_table <- function(data_gene_drugs_table = gene_drugs_table,
                                  input = list()) {
  make_gene_drugs_table_raw <- function() {
    unnested_table <-
      data_gene_drugs_table %>%
      dplyr::filter_all(dplyr::any_vars(fav_gene %in% input$content)) %>% #this ensures no error out if combo of gene found + not found
      tidyr::unnest(data) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(fav_drug)
    return(unnested_table)
  }
  #error handling
  tryCatch(make_gene_drugs_table_raw(),
           error = function(e){
             #make empty table equivalent to returned table
             return(data_gene_drugs_table %>%
                      tidyr::unnest(cols = c(data)) %>%
                      dplyr::slice(0))
           })
}

#' Cell Drugs Table
#'
#' @importFrom magrittr %>%
#'
#' @export
make_cell_drugs_table <- function(data_universal_prism_long = universal_prism_long,
                                  data_cell_expression_meta = cell_expression_meta,
                                  input = list()) {
  make_cell_drugs_table_raw <- function() {
    cell_table <-
      data_universal_prism_long %>%
      dplyr::rename(X1 = 1) %>%
      dplyr::left_join(data_cell_expression_meta, by = "X1") %>%
      dplyr::select(cell_line, name, log2fc) %>%
      dplyr::filter(cell_line %in% input$content) %>%
      dplyr::arrange(dplyr::desc(abs(log2fc)))
    return(cell_table)
  }
  #error handling
  tryCatch(make_cell_drugs_table_raw(),
           error = function(e){
             return(data_universal_prism_long %>%
                      dplyr::rename(X1 = 1) %>%
                      dplyr::left_join(data_cell_expression_meta, by = "X1") %>%
                      dplyr::select(cell_line, name, log2fc) %>%
                      dplyr::slice(0)
             )
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
#' make_metabolite_table(input = list(type = "gene", content = "ROCK2"))
#' make_metabolite_table(input = list(type = "gene", content = c("ROCK1", "ROCK2")))
#' make_metabolite_table(input = list(type = "compound", content = "NAD"))
#' make_metabolite_table(input = list(type = "cell", content = "HEPG2"))
#' \dontrun{
#' make_metabolite_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_metabolite_table <- function(data_compound_hmdb_proteins = compound_hmdb_proteins,
                                  data_compound_hmdb_metabolites = compound_hmdb_metabolites,
                                  data_cell_metabolites = cell_metabolites,
                                  data_cell_expression_names = cell_expression_names,
                                  input = list()) {
  make_metabolite_table_raw <- function() {
    if(input$type == "gene") {
      table_complete <-
        data_compound_hmdb_proteins %>% #gene centric df
        dplyr::filter_all(dplyr::any_vars(fav_gene %in% input$content)) %>%
        tidyr::unnest(data_original) %>%
        dplyr::ungroup() %>%
        dplyr::select(gene_name, metabolite_name, gene_accession, metabolite_accession) %>%
        dplyr::arrange(gene_name, metabolite_name)
    } else if(input$type == "compound") {
      table_complete <-
        data_compound_hmdb_metabolites %>% #metabolite centric df
        dplyr::filter_all(dplyr::any_vars(fav_metabolite %in% input$content)) %>%
        tidyr::unnest(data) %>%
        dplyr::ungroup() %>%
        dplyr::select(metabolite_name, gene_name, gene_accession, metabolite_accession) %>%
        dplyr::arrange(metabolite_name, gene_name)
    } else if(input$type == "cell") {
      table_complete <-
        data_cell_metabolites %>%
        dplyr::left_join(data_cell_expression_names, by = c("DepMap_ID" = "X1")) %>%
        dplyr::filter(cell_line %in% input$content) %>%
        dplyr::select(-CCLE_ID, -DepMap_ID, -lineage, -lineage_subtype) %>%
        tidyr::pivot_longer(cols = !cell_line, names_to = "metabolite") %>%
        dplyr::arrange(-value)
    } else {
      stop("declare your type") }
    return(table_complete)
  }
  #error handling
  tryCatch(make_metabolite_table_raw(),
           error = function(e){
             return(
               data_compound_hmdb_proteins %>%
                 tidyr::unnest(data_original) %>%
                 dplyr::ungroup() %>%
                 dplyr::select(gene_name, metabolite_name, gene_accession, metabolite_accession) %>%
                 dplyr::slice(0))
           })
}
