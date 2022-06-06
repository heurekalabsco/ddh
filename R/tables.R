# HELPER -----
make_empty_table <- function() {
  empty <- as.data.frame("No data available")
  return(empty)
}

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
#' @export
#' @examples
#' make_pathway_list(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' \dontrun{
#' make_pathway_list(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_pathway_list <- function(table_name = pathways,
                              input = list()) { #makes a subtable of pathways that contains your gene query
  if (is.null(input$content)) {
    return (NULL)
  }
  present <- function(list, query){ #is the gene present?
    y <- unlist(list, use.names = FALSE)
    any(str_detect(y, query))
  }
  filtered_table <-
    table_name %>%
    filter(map_lgl(table_name$data, ~present(list = .x, query = input$content)) == TRUE)
  return(filtered_table)
}
#TEST
#make_pathway_list(table_name = pathways, input = list(content = "FXR1"))

# pathway genes makes a gene table of a pathway query
#' Pathway Genes Table
#'
#' \code{make_pathway_genes} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a pathway genes table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a pathway genes table. If an error is thrown, then will return an empty table.
#'
#' @export
#' @examples
#' make_pathway_genes(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' \dontrun{
#' make_pathway_genes(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_pathway_genes <- function(table_name = pathways, table_join = gene_summary, go_id) {
  pathway_table <-
    table_name %>%
    dplyr::filter(go %in% go_id) %>%
    tidyr::unnest(data) %>%
    dplyr::left_join(gene_summary, by = c("gene" = "approved_symbol")) %>%
    dplyr::select(gene, approved_name, aka)
  return(pathway_table)
}

##compounds
#' Compound Table
#'
#' \code{make_compound_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a compound Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a compound Table. If an error is thrown, then will return an empty table.
#'
#' @export
#' @examples
#' make_compound_table(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' \dontrun{
#' make_compound_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_compound_table <- function(data_table = prism_cor_nest,
                                join_table = prism_names,
                                input = list(),
                                top = TRUE) {
  make_compound_table_raw <- function() {

    upper <- prism_cor_upper
    lower <- prism_cor_lower

    table_complete <-
      data_table %>%
      dplyr::filter_all(any_vars(fav_drug %in% input$content)) %>%
      tidyr::unnest("data") %>%
      dplyr::ungroup() %>%
      {if (top == TRUE) filter(., r2 > upper) else filter(., r2 < lower)} %>% #mean +/- 3sd
      dplyr::arrange(desc(r2)) %>%
      dplyr::left_join(join_table, by = "name") %>%
      dplyr::select(1:4)
    return(table_complete)
  }
  #error handling
  tryCatch(make_compound_table_raw(),
           error = function(x){make_empty_table()})
}

#test
#make_compound_table(input = list(compound = "aspirin"), top = TRUE)
#make_compound_table(input = list(compound = "aspirin"), top = FALSE)

# PUBMED TABLE -----
#' Pubmed Table
#'
#' \code{make_pubmed_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a pubmed Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a pubmed Table. If an error is thrown, then will return an empty table.
#'
#' @export
#' @examples
#' make_pubmed_table(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' \dontrun{
#' make_pubmed_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_pubmed_table <- function(pubmed_data = pubmed,
                              input = list()) {
  make_pubmed_table_raw <- function() {
    if(input$type == "gene") {
      pubmed_table <-
        pubmed_data %>%
        dplyr::filter_all(any_vars(name %in% input$content)) %>%
        tidyr::unnest(data) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(as.numeric(pmid))
    } else if(input$type == "compound") {
      pubmed_table <-
        pubmed_data %>%
        dplyr::filter_all(any_vars(name %in% input$content)) %>%
        tidyr::unnest(data) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(as.numeric(pmid))
    } else if(input$type == "cell" | input$type == "lineage" | input$type == "lineage_subtype" | input$type == "cell_list") {
      pubmed_table <-
        pubmed_data %>%
        dplyr::filter(name %in% input$content) %>%
        tidyr::unnest(data) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(as.numeric(pmid))
    }else {
      stop("declare your type") }
    return(pubmed_table)
  }
  #error handling
  tryCatch(make_pubmed_table_raw(),
           error = function(x){make_empty_table()})
}

#make_pubmed_table(input = list(type = "gene", content = "DXTL3"))
#make_pubmed_table(input = list(type = "gene", content = c("ROCK1", "ROCK2"))) %>% count()
#make_pubmed_table(input = list(type = "gene", content = "ROCK"))

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
#' @export
#' @examples
#' make_cellanatogram_table(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' \dontrun{
#' make_cellanatogram_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_cellanatogram_table <- function(cellanatogram_data = subcell,
                                     input = list()) { #change to input=list()
  make_cellanatogram_table_raw <- function() {
    cellanatogram_data %>%
      dplyr::filter_all(any_vars(gene_name %in% input$content)) %>%
      dplyr::filter(!is.na(type)) %>%
      dplyr::add_count(main_location) %>%
      dplyr::transmute(Gene = gene_name, Reliability = reliability, Location = main_location, Count = as_factor(n)) %>%
      dplyr::arrange(desc(Count))
  }
  #error handling
  tryCatch(make_cellanatogram_table_raw(),
           error = function(x){make_empty_table()})
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
#' @export
#' @examples
#' make_expression_table(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' \dontrun{
#' make_expression_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_expression_table <- function(expression_data = expression_long,
                                  expression_join = expression_names,
                                  input = list(),
                                  var = "gene") { #you are so slow
  make_expression_table_raw <- function() {
    if (var == "gene") {
      table_data <-
        expression_data %>%
        dplyr::select(any_of(c("X1", "gene", "gene_expression"))) %>%
        dplyr::rename("expression_var" = "gene_expression")
    } else if (var == "protein") {
      table_data <-
        expression_data %>%
        dplyr::select(any_of(c("X1", "gene", "protein_expression"))) %>%
        dplyr::rename("expression_var" = "protein_expression")
    } else {
      stop("delcare your variable")
    }
    if (input$type == "gene") {
      table_data %>%
        dplyr::filter_all(any_vars(gene %in% input$content)) %>%
        dplyr::arrange(desc(.[3])) %>%
        tidyr::pivot_wider(names_from = gene, values_from = expression_var) %>%
        dplyr::left_join(expression_join, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::select(cell_line, lineage, lineage_subtype, everything()) %>%
        dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
        dplyr::rename("Cell Line" = "cell_line", "Lineage" = "lineage", "Subtype" = "lineage_subtype")
    } else if (input$type == "cell") {
      table_data %>%
        dplyr::arrange(desc(.[3])) %>%
        dplyr::left_join(expression_join, by = "X1") %>%
        dplyr::select(-X1, -lineage, -lineage_subtype) %>%
        dplyr::filter_all(any_vars(cell_line %in% input$content)) %>%
        tidyr::pivot_wider(names_from = cell_line, values_from = expression_var) %>%
        dplyr::select(gene, everything()) %>%
        dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
        dplyr::rename("Gene" = "gene")
    }
  }
  #error handling
  tryCatch(make_expression_table_raw(),
           error = function(x){make_empty_table()})
}

#make_expression_table(input = list(type = "gene", subtype = "gene", content = "ROCK1"))
#make_expression_table(input = list(type = "gene", subtype = "gene", content = "HSP90B1"))

# HUMAN ANATOGRAM TABLES -----
#' Humananatogram Table
#'
#' \code{make_humananatogram_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a humananatogram Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a humananatogram Table. If an error is thrown, then will return an empty table.
#'
#' @export
#' @examples
#' make_humananatogram_table(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' \dontrun{
#' make_humananatogram_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_humananatogram_table <- function(humananatogram_data = tissue,
                                      input = list()) {
  make_humananatogram_table_raw <- function() {
    humananatogram_data %>%
      dplyr::filter_all(any_vars(gene_name %in% input$content)) %>%
      dplyr::filter(!is.na(value)) %>%
      dplyr::mutate(organ = stringr::str_replace_all(organ, "_", " "),
                    organ = stringr::str_to_title(organ)) %>%
      dplyr::arrange(desc(value)) %>%
      dplyr::select(-gene) %>%
      tidyr::pivot_wider(names_from = gene_name, values_from = value, names_sort = TRUE) %>%
      dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
      dplyr::rename("Organ" = "organ")
  }
  #error handling
  tryCatch(make_humananatogram_table_raw(),
           error = function(x){make_empty_table()})
}

#make_humananatogram_table(input = list(type = "gene", subtype = "gene", content = "ROCK1"))

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
#' @export
#' @examples
#' make_clustering_table(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' \dontrun{
#' make_clustering_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_clustering_table <- function(cluster_data = sequence_clusters,
                                  signature_data = signatures,
                                  cluster_names = protein_cluster_names,
                                  input = list(),
                                  cluster = FALSE,
                                  show_signature = FALSE) {

  make_clustering_table_raw <- function() {

    # sequence_data_clean <- cluster_data %>%
    #   dplyr::filter(clust != 0)

    # sequence_data_clean <- cluster_data

    query_clust <- cluster_data %>%
      dplyr::filter(gene_name %in% input$content) %>%
      dplyr::pull(clust) %>%
      unique()

    cluster_table <- cluster_data %>%
      dplyr::left_join(signature_data, by = "uniprot_id") %>%
      dplyr::select(uniprot_id, gene_name.x, protein_name,
                    clust) %>%
      dplyr::filter(clust %in% query_clust) %>%
      dplyr::rename(gene_name = 2) %>%
      dplyr::left_join(cluster_names, by = "clust") %>%
      dplyr::relocate(cluster_name, .after = clust)

    if(!cluster) {
      cluster_table <- cluster_table %>%
        dplyr::filter(gene_name %in% input$content)
    }

    if(show_signature) {
      cluster_table <- cluster_table %>%
        dplyr::left_join(signature_data %>%
                           dplyr::select(uniprot_id, A:U) %>%
                           mutate_at(vars(A:U), ~ round(., 2)),
                         by = "uniprot_id")
    }

    if(nrow(cluster_table) == 0) {
      stop("Unable to cluster this protein by its amino acid sequence so far...")
    }

    return(cluster_table)

  }

  #error handling
  tryCatch(make_clustering_table_raw(),
           error = function(x){make_empty_table()})
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
#' @export
#' @examples
#' make_clustering_enrichment_table(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' \dontrun{
#' make_clustering_enrichment_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_clustering_enrichment_table <- function(cluster_data = sequence_clusters,
                                             input = list(),
                                             enrichment_tables = enriched_clusters,
                                             ontology = "BP") {

  make_clustering_enrichment_table_raw <- function() {

    sequence_data_clean <-
      cluster_data %>%
      dplyr::filter(clust != 0)

    query_clust <-
      sequence_data_clean %>%
      dplyr::filter(gene_name %in% input$content) %>%
      dplyr::pull(clust) %>%
      unique()

    if(length(query_clust) != 1) {
      stop("Select only one cluster")
    }

    enrichment_table <-
      enrichment_tables %>%
      bind_rows() %>%
      dplyr::filter(cluster == query_clust) %>%
      dplyr::filter(ont == ontology) %>%
      dplyr::select(-ont, -cluster)

    if(nrow(enrichment_table) == 0) {
      stop("No enriched terms for this cluster and ontology...")
    }

    return(enrichment_table)
  }

  #error handling
  tryCatch(make_clustering_enrichment_table_raw(),
           error = function(x){make_empty_table()})
}

# DEPENDENCY TABLES -----
#' Dep Table
#'
#' \code{make_dep_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a dep Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a dep Table. If an error is thrown, then will return an empty table.
#'
#' @export
#' @examples
#' make_dep_table(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' \dontrun{
#' make_dep_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_dep_table <- function(achilles_data = achilles_long,
                           prism_data = prism_long,
                           prism_summary = prism_names,
                           expression_data = expression_names,
                           input = list(),
                           var = "gene") { #Variable to determine type of dependency table for cell lines
  make_dep_table_raw <- function() {
    if(input$type == "gene") {
      table_data <-
        achilles_data %>%
        dplyr::filter_all(any_vars(gene %in% input$content)) %>%
        tidyr::pivot_wider(names_from = gene, values_from = dep_score) %>%
        dplyr::left_join(expression_data, by = "X1") %>%
        dplyr::select(cell_line, lineage, contains(input$content)) %>%
        dplyr::rename("Cell Line" = "cell_line", "Lineage" = "lineage")
    } else if(input$type == "compound") {
      table_data <-
        prism_data %>%
        dplyr::filter_all(any_vars(name %in% input$content)) %>%
        tidyr::pivot_wider(names_from = name, values_from = log2fc) %>%
        dplyr::left_join(expression_data, by = c("x1" = "X1")) %>%
        dplyr::select(cell_line, lineage, contains(input$content)) %>%
        dplyr::rename("Cell Line" = "cell_line", "Lineage" = "lineage")
    } else {
      if(var == "gene"){
        table_data <-
          achilles_data %>%
          dplyr::left_join(expression_data, by = "X1") %>%
          dplyr::filter(cell_line %in% input$content) %>%
          tidyr::pivot_wider(names_from = cell_line, values_from = dep_score) %>%
          dplyr::rename("dep_score" = as.character(input$content)) %>%
          dplyr::left_join(gene_summary, by = c("gene" = "approved_symbol")) %>%
          dplyr::select(gene, approved_name, dep_score, contains(input$content)) %>%
          dplyr::mutate(unique_essential = gene %in% unique_essential_genes$gene, common_essential = gene %in% common_essentials$gene)
      }
      else if (var == "drug"){
        table_data <-
          prism_data %>%
          dplyr::left_join(expression_data, by = c("x1" = "X1")) %>%
          dplyr::left_join(prism_summary, by = "name") %>%
          dplyr::filter(cell_line %in% input$content) %>%
          tidyr::pivot_wider(names_from = cell_line, values_from = log2fc) %>%
          dplyr::rename("log2fc" = as.character(input$content)) %>%
          dplyr::select(name, moa, log2fc, contains(input$content)) %>%
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
           error = function(x){make_empty_table()})
}

#test
#make_dep_table(input = list(type = "gene", query = "ROCK1", content = "ROCK1"))
#make_dep_table(input = list(type = "gene", query = "ROCK1", content = c("ROCK1", "ROCK2")))
#make_dep_table(input = list(type = "compound", query = "aspirin", content = "aspirin"))
#make_dep_table(input = list(type = "cell", query = "HEPG2", content = "HEPG2"))

##co-essentiality-----
##similar
#' Top Table
#'
#' \code{make_top_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a top Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a top Table. If an error is thrown, then will return an empty table.
#'
#' @export
#' @examples
#' make_top_table(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' \dontrun{
#' make_top_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_top_table <- function(toptable_data = master_top_table,
                           toptable_data_cell = master_top_table_cell_line,
                           upper = achilles_upper,
                           gls = FALSE,
                           input = list()) {
  make_top_table_raw <- function() {
    if(input$type == "gene") {
      table_data <-
        toptable_data %>%
        dplyr::filter_all(any_vars(fav_gene %in% input$content))
    } else if(input$type == "cell") {
      table_data <-
        toptable_data_cell %>%
        dplyr::filter_all(any_vars(fav_cell %in% input$content))
    } else {
      stop("delcare your type")
    }
    if(nrow(table_data) != 0) {
      if(gls){
        table_complete <-
          table_data %>%
          tidyr::unnest(data) %>%
          filter(GLSpvalue < 0.05 & !duplicated(gene)) %>%
          dplyr::arrange(-desc(GLSpvalue))
      } else {
        table_complete <-
          table_data %>%
          tidyr::unnest(data) %>%
          filter(r2 > upper & !duplicated(gene)) %>%
          dplyr::arrange(desc(r2))
      }
    } else {
      table_complete <-
        table_data
    }
    return(table_complete)
  }
  #error handling
  tryCatch(make_top_table_raw(),
           error = function(x){make_empty_table()})
}
#test
#make_top_table(input = list(type = "gene", content = "ROCK1"))
#make_top_table(input = list(type = "gene", content = "ROCK"))
#make_top_table(input = list(type = "cell", content = "MCC13"))

##dissimilar
#' Bottom Table
#'
#' \code{make_bottom_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a bottom Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a bottom Table. If an error is thrown, then will return an empty table.
#'
#' @export
#' @examples
#' make_bottom_table(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' \dontrun{
#' make_bottom_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_bottom_table <- function(bottomtable_data = master_bottom_table,
                              bottomtable_data_cell = master_bottom_table_cell_line,
                              lower = achilles_lower,
                              gls = FALSE,
                              input = list()) {
  make_bottom_table_raw <- function() {
    if(input$type == "gene") {
      table_data <-
        bottomtable_data %>%
        dplyr::filter_all(any_vars(fav_gene %in% input$content))
    } else if(input$type == "cell") {
      table_data <-
        bottomtable_data_cell %>%
        dplyr::filter_all(any_vars(fav_cell %in% input$content))
    } else {
      stop("delcare your type")
    }
    if(nrow(table_data) != 0) {
      if(gls){
        table_complete <-
          table_data %>%
          tidyr::unnest(data) %>%
          filter(GLSpvalue < 0.05 & !duplicated(gene)) %>%
          dplyr::arrange(-desc(GLSpvalue))
      } else {
        table_complete <-
          table_data %>%
          tidyr::unnest(data) %>%
          filter(r2 < lower & !duplicated(gene)) %>%
          dplyr::arrange(desc(r2))
      }
    } else {
      table_complete <-
        table_data
    }
    return(table_complete)
  }
  #error handling
  tryCatch(make_bottom_table_raw(),
           error = function(x){make_empty_table()})
}

#test
#make_bottom_table(input = list(type = "gene", content = "ROCK1"))
#make_bottom_table(input = list(type = "gene", content = "QARS1"))

##censor-----
#censor is used to remove genes from similarity table that are garbage (too many associations)
censor <- function(top_table,
                   censor_data = censor_genes,
                   choice = FALSE,
                   greater_than) {
  if(choice == TRUE){
    censor_data <-
      censor_data %>%
      dplyr::filter(num_sim > greater_than) #get the genes that have too many associations

    censored_table <-
      top_table %>%
      dplyr::filter(!gene %in% censor_data$genes) #filter NOT in (out) the too high list

    return(censored_table)
  }
  return(top_table)
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
#' @export
#' @examples
#' make_enrichment_top(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' \dontrun{
#' make_enrichment_top(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_enrichment_top <- function(enrichmenttop_data = master_positive,
                                input = list()) {
  make_enrichment_top_raw <- function() {
    enrichmenttop_data %>%
      dplyr::filter_all(any_vars(fav_gene %in% input$content)) %>%
      tidyr::unnest(data) %>%
      dplyr::select(fav_gene, enrichr, Term, Overlap, Adjusted.P.value, Combined.Score, Genes) %>%
      dplyr::arrange(Adjusted.P.value) %>%
      dplyr::rename("Query" = "fav_gene", "Gene Set" = "enrichr", "Gene List" = "Term", "Adjusted p-value" = "Adjusted.P.value", "Combined Score" = "Combined.Score") #"Overlap", "Genes"
  }
  #error handling
  tryCatch(make_enrichment_top_raw(),
           error = function(x){make_empty_table()})
}

# make_enrichment_top(input = list(content = c("ROCK1")))
# make_enrichment_top(input = list(type = "gene", query = "ROCK1", content = c("ROCK1")))

#' Enrichment Bottom Table
#'
#' \code{make_enrichment_bottom} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a enrichment bottom table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a enrichment bottom table. If an error is thrown, then will return an empty table.
#'
#' @export
#' @examples
#' make_enrichment_bottom(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' \dontrun{
#' make_enrichment_bottom(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_enrichment_bottom <- function(enrichmentbottom_data = master_negative,
                                   input = list()) {
  make_enrichment_bottom_raw <- function() {
    enrichmentbottom_data %>%
      dplyr::filter_all(any_vars(fav_gene %in% input$content)) %>%
      tidyr::unnest(data) %>%
      dplyr::select(fav_gene, enrichr, Term, Overlap, Adjusted.P.value, Combined.Score, Genes) %>%
      dplyr::arrange(Adjusted.P.value) %>%
      dplyr::rename("Query" = "fav_gene", "Gene Set" = "enrichr", "Gene List" = "Term", "Adjusted p-value" = "Adjusted.P.value", "Combined Score" = "Combined.Score") #"Overlap", "Genes"
  }
  #error handling
  tryCatch(make_enrichment_bottom_raw(),
           error = function(x){make_empty_table()})
}

# make_enrichment_bottom(input = list(content = c("ROCK1")))

# CELL SIMILARITY TABLE -----
#' Cell Sim Table
#'
#' \code{make_cell_sim_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a cell sim Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a cell sim Table. If an error is thrown, then will return an empty table.
#'
#' @export
#' @examples
#' make_cell_sim_table(input = list(type = "cell", content = "HEL"))
#' make_cell_sim_table(input = list(type = "cell", content = c("HEL", "LS513")))
#' \dontrun{
#' make_cell_sim_table(input = list(type = "cell", content = "HEL"))
#' }
make_cell_sim_table <- function(cell_sims = cell_line_combs_gam_p,
                                cell_metadata = expression_names,
                                sim_type = "bio",
                                pval_ct = 0.01,
                                fdr_ct = 0.05,
                                input = list()) {

  make_cell_sim_table_raw <- function() {
    cell_sim_table <-
      cell_sims %>%
      dplyr::filter(cell1_name %in% input$content | cell2_name %in% input$content)

    # Swap cols
    cell_sim_table[cell_sim_table$cell2_name %in% input$content, c("cell1_name", "cell2_name")] <-
      cell_sim_table[cell_sim_table$cell2_name %in% input$content, c("cell2_name", "cell1_name")]

    # Linear and non-linear associations are provided in the same table.
    # However, p-values should be interpreted separately.
    # Important to highlight in "methods".
    if(sim_type == "bio") {
      cell_sim_table <- cell_sim_table %>%
        dplyr::mutate(pvalue = ifelse(type_bio == "LM", lm_p, gam_p),
                      fdr = ifelse(type_bio == "LM", lm_fdr, gam_fdr),
                      association_type = ifelse(type_bio == "LM", "Linear", "Non-linear")) %>%
        dplyr::select(cell1_name, cell2_name, pvalue, fdr, association_type) %>%
        dplyr::filter(pvalue < pval_ct & fdr < fdr_ct) %>%
        dplyr::rename(cell1 = cell1_name, cell2 = cell2_name)

    } else {
      cell_sim_table <- cell_sim_table %>%
        dplyr::mutate(pvalue = ifelse(type_stat == "LM", lm_p, gam_p),
                      fdr = ifelse(type_stat == "LM", lm_fdr, gam_fdr),
                      association_type = ifelse(type_stat == "LM", "Linear", "Non-linear")) %>%
        dplyr::select(cell1_name, cell2_name, pvalue, fdr, association_type) %>%
        dplyr::filter(pvalue < pval_ct & fdr < fdr_ct) %>%
        dplyr::rename(cell1 = cell1_name, cell2 = cell2_name)

    }

    cell_sim_table <- cell_sim_table %>%
      # dplyr::left_join(cell_metadata %>%
      #                    dplyr::select(cell_line, lineage, lineage_subtype),
      #                  by = c("cell1" = "cell_line")
      #                  ) %>%
      # dplyr::rename(lineage_cell1 = lineage,
      #               lineage_subtype_cell1 = lineage_subtype) %>%
      dplyr::left_join(cell_metadata %>%
                         dplyr::select(cell_line, lineage, lineage_subtype),
                       by = c("cell2" = "cell_line")
      ) %>%
      dplyr::rename(lineage_cell2 = lineage,
                    lineage_subtype_cell2 = lineage_subtype) %>%
      dplyr::select(cell1, cell2,
                    # lineage_cell1, lineage_subtype_cell1,
                    lineage_cell2, lineage_subtype_cell2,
                    pvalue, fdr,
                    association_type)

    return(cell_sim_table)
  }

  # error handling
  tryCatch(make_cell_sim_table_raw(),
           error = function(x){make_empty_table()})
}

# DRUG TABLES -----
#' Drug Genes Cor Table
#'
#' \code{make_drug_genes_cor_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a drug genes cor Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a drug genes cor Table. If an error is thrown, then will return an empty table.
#'
#' @export
#' @examples
#' make_drug_genes_table(input = list(content = "aspirin"))
#' \dontrun{
#' make_drug_genes_table(input = list(content = "aspirin"))
#' }
make_drug_genes_cor_table <- function(table_data = drug_genes_cor_table, drug) {
  make_drug_genes_cor_table_raw <- function() {
    unnested_table <-
      table_data %>%
      dplyr::filter_all(any_vars(fav_drug %in% drug)) %>%
      tidyr::unnest(data) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(desc(r2))
    return(unnested_table)
  }
  tryCatch(make_drug_genes_cor_table_raw(),
           error = function(x){make_empty_table()})
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
#' @export
#' @examples
#' make_gene_drugs_cor_table(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' \dontrun{
#' make_gene_drugs_cor_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_gene_drugs_cor_table <- function(table_data = gene_drugs_cor_table,
                                      input = list()) {
  make_gene_drugs_cor_table_raw <- function() {
    unnested_table <-
      table_data %>%
      dplyr::filter_all(any_vars(fav_gene %in% input$content)) %>%
      tidyr::unnest(data) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(desc(r2))
    return(unnested_table)
  }
  tryCatch(make_gene_drugs_cor_table_raw(),
           error = function(x){make_empty_table()})
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
#' @export
#' @examples
#' make_drug_genes_table(drug = "aspirin")
#' make_drug_genes_table(drug = "ibuprofen")
#' make_drug_genes_table(drug = c("aspirin", "ibuprofen"))
#' \dontrun{
#' make_drug_genes_table(drug = "aspirin")
#' }
make_drug_genes_table <- function(table_data = drug_genes_table,
                                  drug) {
  make_drug_genes_table_raw <- function() {
    unnested_table <-
      table_data %>%
      dplyr::filter_all(any_vars(fav_drug %in% drug)) %>%
      tidyr::unnest(data) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(fav_gene)
    return(unnested_table)
  }
  #error handling
  tryCatch(make_drug_genes_table_raw(),
           error = function(x){make_empty_table()})
}

#' Gene Drugs Table
#'
#' \code{make_gene_drugs_table} returns an image of ...
#'
#' This is a table function that takes a gene name and returns a gene drugs Table
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a gene drugs Table. If an error is thrown, then will return an empty table.
#'
#' @export
#' @examples
#' make_gene_drugs_table(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_gene_drugs_table_raw(input = list(content = c("ROCK1", "RDX")))
#' make_gene_drugs_table_raw(input = list(content = c("ROCK1", "RDX", "ROCK")))
#' \dontrun{
#' make_gene_drugs_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_gene_drugs_table <- function(table_data = gene_drugs_table,
                                  input = list()) {
  make_gene_drugs_table_raw <- function() {
    unnested_table <-
      table_data %>%
      dplyr::filter_all(any_vars(fav_gene %in% input$content)) %>% #this ensures no error out if combo of gene found + not found
      tidyr::unnest(data) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(fav_drug)
    return(unnested_table)
  }
  #error handling
  tryCatch(make_gene_drugs_table_raw(),
           error = function(x){make_empty_table()})
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
#' @export
#' @examples
#' make_metabolite_table(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_metabolite_table(input = list(type = "gene", content = "ROCK2"))
#' make_metabolite_table(input = list(type = "gene", content = c("ROCK1", "ROCK2")))
#' make_metabolite_table(input = list(type = "compound", content = "aspirin")
#' make_metabolite_table(input = list(type = "cell", query = "HEPG2", content = "HEPG2"))
#' \dontrun{
#' make_metabolite_table(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_metabolite_table <- function(protein_data = hmdb_proteins,
                                  metabolite_data = hmdb_metabolites,
                                  input = list()) {
  make_metabolite_table_raw <- function() {
    if(input$type == "gene") {
      table_complete <-
        protein_data %>% #gene centric df
        dplyr::filter_all(any_vars(fav_gene %in% input$content)) %>%
        tidyr::unnest(data_original) %>%
        dplyr::ungroup() %>%
        dplyr::select(gene_name, metabolite_name, gene_accession, metabolite_accession) %>%
        dplyr::arrange(gene_name, metabolite_name)
    } else if(input$type == "compound") {
      table_complete <-
        metabolite_data %>% #metabolite centric df
        dplyr::filter_all(any_vars(fav_metabolite %in% input$content)) %>%
        tidyr::unnest(data) %>%
        dplyr::ungroup() %>%
        dplyr::select(metabolite_name, gene_name, gene_accession, metabolite_accession) %>%
        dplyr::arrange(metabolite_name, gene_name)
    } else if(input$type == "cell") {
      table_complete <-
        metabolite_data %>%
        dplyr::left_join(expression_names, by = c("DepMap_ID" = "X1")) %>%
        dplyr::filter(cell_line %in% input$content) %>%
        dplyr::select(-CCLE_ID, -DepMap_ID, -lineage, -lineage_subtype) %>%
        pivot_longer(cols = !cell_line, names_to = "metabolite") %>%
        arrange(-value)
    } else {
      stop("declare your type") }
    return(table_complete)
  }
  #error handling
  tryCatch(make_metabolite_table_raw(),
           error = function(x){make_empty_table()})
}
