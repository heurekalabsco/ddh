#' Make Gene Summary
#'
#' The make_gene_summary function takes a gene as an input and returns summary text.
#'
#' @param input A list containing a content variable.
#' @param var Variable that determines which text is returned
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_summary_gene(input = list(content = "ROCK1"), var = "id")
#' make_summary_gene(input = list(content = "ROCK1"), var = "name")
#' make_summary_gene(input = list(content = "ROCK1"), var = "summary")
#' make_summary_gene(input = list(content = c("ROCK1", "ROCK2")), var = "name")
make_summary_gene <- function(input = list(),
                              var = "id") {
  if (is.null(input$content)) {
    return (NULL)
  }
  gene_summary_var <-
    get_data_object(object_name = input$content,
                    dataset_name = "universal_gene_summary",
                    pivotwider = TRUE) %>%
    dplyr::mutate(across(contains(c("count", "rank")), as.numeric)) %>%
    dplyr::rename(name = approved_name, summary = entrez_summary) %>%
    dplyr::pull(var) #any column name

  return(gene_summary_var)
}

#' Make Pathway Summary
#'
#' The make_pathway_summary function takes an id as a character input and returns summary information of the pathway.
#'
#' @param input A list containing a content variable.
#' @param var Variable that determines which text is returned
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_summary_pathway(input = list(query = 5887), var = "gs_name")
#' make_summary_pathway(input = list(query = 5887), var = "gs_cat")
#' make_summary_pathway(input = list(query = 5887), var = "gs_description")
#' make_summary_pathway(input = list(query = 5887), var = "pathway_size")
make_summary_pathway <- function(input = list(),
                                 var = "gs_description") {
  if (is.null(input$query)) {
    return (NULL)
  }
  pathway_summary_var <-
    get_content("universal_pathways", dataset = TRUE) %>%
    dplyr::filter(gs_id %in% input$query) %>%
    dplyr::mutate(gs_name = sub("^[^_]*_", "", gs_name)) %>%
    dplyr::mutate(gs_name = gsub("_", " ", gs_name),
                  gs_name = stringr::str_to_title(gs_name)) %>%
    dplyr::select(dplyr::all_of(var)) %>%
    dplyr::pull(1) %>% #any column name
    unique()
  return(pathway_summary_var)
}

#' Make Summary Text
#'
#' This function takes one or many genes as inputs and returns summary list key summary information. This is a solution for make_summary_gene to make summaries of more than one gene.
#'
#' @param input A list containing a content variable.
#' @param var Variable that determines which text is returned
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_summary_text(input = list(type = "gene", content = c("ROCK1", "ROCK2")))
#' make_summary_text(input = list(type = "cell", content = c("HEL", "HEPG2")))
#' make_summary_text(input = list(type = "gene", subtype = "pathway", query = 5887))
make_summary_text <- function(input = list(),
                              summary_len = 40, # number of WORDS
                              ...) {
  if (is.null(input$query) & is.null(input$content)) {
    return (NULL)
  }

  if (input$type == "gene") {
    if (input$subtype == "pathway") {
      pathway_summary_var <-
        get_content("universal_pathways", dataset = TRUE) %>%
        dplyr::filter(gs_id %in% input$query) %>%
        dplyr::mutate(gs_name = sub("^[^_]*_", "", gs_name)) %>%
        dplyr::mutate(gs_name = gsub("_", " ", gs_name),
                      gs_name = stringr::str_to_title(gs_name))

      valid_summaries <- glue::glue("<div><h3>Pathway: {pathway_summary_var$gs_name[1]}</h3></div>
                                    <div><b>Gene Set ID: </b>{pathway_summary_var$gs_id[1]}</div>
                                    <div><b>Category: </b>{pathway_summary_var$gs_cat[1]}</div>
                                    <div><b>Subategory: </b>{pathway_summary_var$gs_subcat[1]}</div>
                                    <div><b>Size: </b>{pathway_summary_var$pathway_size[1]}</div>
                                    <div><b>URL: </b><a href='{pathway_summary_var$gs_url[1]}' target='_blank'>{pathway_summary_var$gs_url[1]}</a></div>
                                    <div><b>Genes: </b></div>
                                    <div><p>{paste0(pathway_summary_var$human_gene_symbol, collapse = ', ')}</p></div>
                                    <div><b>Description: </b></div>
                                    <div><p>{pathway_summary_var$gs_description[1]}</p></div>
                                    ") %>%
        htmltools::HTML()

    } else {
      data_gene_location <-
        get_data_object(object_names = input$content,
                        dataset_name = "gene_location",
                        pivotwider = TRUE)

      custom_list <-
        get_data_object(object_name = input$content,
                        dataset_name = "universal_gene_summary",
                        pivotwider = TRUE) %>%
        dplyr::left_join(data_gene_location, by = "id") %>%
        dplyr::mutate(across(contains(c("count", "rank")), as.numeric))

      custom_list[custom_list == "NA"] <- NA
      custom_list[custom_list == ""] <- NA
      custom_list <-
        custom_list %>%
        dplyr::mutate_all(~ ifelse(is.na(.), "No info.", .))

      if (length(input$content) == 1) {
        valid_summaries <- glue::glue("<div><h3>{custom_list$id}: {custom_list$approved_name}</h3></div>
                                    <div><b>Entrez ID: </b><a href='https://www.ncbi.nlm.nih.gov/gene/?term={custom_list$ncbi_gene_id}' target='_blank'>{custom_list$ncbi_gene_id}</a></div>
                                    <div><b>ENSEMBL ID: </b>{custom_list$ensembl_gene_id}</div>
                                    <div><b>Chromosome: </b>{custom_list$chromosome}</div>
                                    <div><b>Coding Sequence Length: </b>{custom_list$cds_length} bp</div>
                                    <div><b>Aka: </b>{custom_list$aka}</div>
                                    <div><b>Description: </b></div>
                                    <div><p>{custom_list$entrez_summary}</p></div>
                                    ") %>%
          htmltools::HTML()
      } else {
        custom_list <-
          custom_list %>%
          dplyr::mutate(entrez_summary = ifelse(stringr::str_count(entrez_summary) >= summary_len,
                                                paste0(gsub(paste0("^((\\w+\\W+){", summary_len, "}\\w+).*$"), "\\1", entrez_summary), " ..."),
                                                entrez_summary)
          )

        custom_list_split <- split(custom_list, custom_list$id)

        tab_fun_gene <- function(custom_list_split){
          summary_tables <- list()
          for (i in names(custom_list_split)){
            tabledata <- custom_list_split[[i]]
            summary_tables[[i]] <- glue::glue("<div><a href='?show=gene&query={tabledata$id}' target='_blank'><h3>{tabledata$id}</a>: {tabledata$approved_name}</h3></div>
                                            <div><b>Entrez ID: </b><a href='https://www.ncbi.nlm.nih.gov/gene/?term={tabledata$ncbi_gene_id}' target='_blank'>{tabledata$ncbi_gene_id}</a></div>
                                            <div><b>ENSEMBL ID: </b>{tabledata$ensembl_gene_id}</div>
                                            <div><b>Description: </b></div>
                                            <div><p>{tabledata$entrez_summary}</p></div>
                                            ")
          }
          return(dplyr::bind_rows(summary_tables) %>%
                   tidyr::unite("text", dplyr::everything(), sep = " ")
          )
        }

        valid_summaries <-
          custom_list_split %>%
          tab_fun_gene() %>%
          dplyr::pull(text) %>%
          htmltools::HTML()
      }
    }
  }

  # } else if (input$type == "cell") {
    # custom_list <-
    #   data_cell_expression_meta %>%
    #   dplyr::filter(cell_line %in% input$content) %>%
    #   dplyr::select(cell_line, lineage, lineage_subtype, age, sex) %>%
    #   dplyr::left_join(data_cell_osaurus %>%
    #                      dplyr::select(name, CC) %>%
    #                      dplyr::rename(cell_line = name),
    #                    by = c("cell_line")) %>%
    #   tidyr::replace_na(list(CC = "NA"))
    #
    # custom_list[custom_list == "NA"] <- NA
    # custom_list[custom_list == ""] <- NA
    # custom_list <- custom_list %>%
    #   dplyr::mutate_all(~ ifelse(is.na(.), "No info.", .))
    #
    # if(length(input$content) == 1) {
    #
    #   valid_summaries <- glue::glue("<div><h3>{custom_list$cell_line} ({custom_list$lineage})</h3></div>
    #                                 <div><b>Lineage: </b>{custom_list$lineage}</div>
    #                                 <div><b>Lineage subtype: </b>{custom_list$lineage_subtype}</div>
    #                                 <div><b>Age: </b>{custom_list$age}</div>
    #                                 <div><b>Sex: </b>{custom_list$sex}</div>
    #                                 <div><b>Description</b></div>
    #                                 <div><p>{custom_list$CC}</p></div>
    #                                 ") %>%
    #     htmltools::HTML()
    # } else {
    #   custom_list <-
    #     custom_list %>%
    #     dplyr::mutate(CC = ifelse(stringr::str_count(CC) >= summary_len,
    #                               paste0(gsub(paste0("^((\\w+\\W+){", summary_len, "}\\w+).*$"), "\\1", CC), " ..."),
    #                               CC)
    #     )
    #
    #   custom_list_split <- split(custom_list, custom_list$cell_line)
    #
    #   tab_fun_cell <- function(custom_list_split){
    #     summary_tables <- list()
    #     for (i in names(custom_list_split)){
    #       tabledata <- custom_list_split[[i]]
    #       summary_tables[[i]] <- glue::glue("<div><a href='?show=cell&query={tabledata$cell_line}' target='_blank'><h3>{tabledata$cell_line} ({tabledata$lineage})</a></h3></div>
    #                                         <div><b>Lineage: </b>{tabledata$lineage}</div>
    #                                         <div><b>Lineage subtype: </b>{tabledata$lineage_subtype}</div>
    #                                         <div><b>Description</b></div>
    #                                         <div><p>{tabledata$CC}</p></div>
    #                                         ")
    #     }
    #     return(dplyr::bind_rows(summary_tables) %>%
    #              tidyr::unite("text", dplyr::everything(), sep = " ")
    #     )
    #   }
    #
    #   valid_summaries <-
    #     custom_list_split %>%
    #     tab_fun_cell() %>%
    #     dplyr::pull(text) %>%
    #     htmltools::HTML()
    # }

  return(valid_summaries)
}

#' Make Protein Summary
#'
#' The make_summary_protein function takes a gene as an input and returns summary text about its protein.
#'
#' @param input A list containing a content variable.
#' @param var Variable that determines which text is returned
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_summary_protein(input = list(type = "gene", content = c("ROCK1")))
#' make_summary_protein(input = list(type = "gene", content = c("ROCK1", "ROCK2")))
make_summary_protein <- function(input = list(),
                                 summary_len = 40, # number of WORDS
                                 ...) {
  if (is.null(input$content)) {
    return (NULL)
  }

  custom_list <-
    get_data_object(object_name = input$content,
                    dataset_name = "universal_proteins",
                    pivotwider = TRUE) %>%
    dplyr::mutate(across(contains(c("mass")), as.numeric))

  custom_list[custom_list == "NA"] <- NA
  custom_list[custom_list == ""] <- NA
  custom_list <- custom_list %>%
    dplyr::mutate_all(~ ifelse(is.na(.), "No info.", .)) %>%
    dplyr::mutate(function_cc = gsub("\\s*\\{[^\\)]+\\}*.", "", function_cc))

  if (length(input$content) == 1) {
    valid_summaries <- glue::glue("<div><h3>{custom_list$id}: {custom_list$protein_name}</h3></div>
                                  <div><b>Uniprot ID: </b><a href='https://www.uniprot.org/uniprot/{custom_list$uniprot_id}' target='_blank'>{custom_list$uniprot_id}</a></div>
                                  <div><b>Enzyme Commission: </b><a href='https://enzyme.expasy.org/EC/{custom_list$ec}' target='_blank'>{custom_list$ec}</a></div>
                                  <div><b>Protein Mass: </b>{custom_list$mass} kDa</div>
                                  <div><b>Description: </b></div>
                                  <div><p>{custom_list$function_cc}</p></div>
                                  ") %>%
      htmltools::HTML()
  } else {
    custom_list <- custom_list %>%
      dplyr::mutate(function_cc = ifelse(stringr::str_count(function_cc) >= summary_len,
                                         paste0(gsub(paste0("^((\\w+\\W+){", summary_len, "}\\w+).*$"), "\\1", function_cc), " ..."),
                                         function_cc)
      )

    custom_list_split <- split(custom_list, custom_list$id)

    tab_fun_gene <- function(custom_list_split) {
      summary_tables <- list()
      for (i in names(custom_list_split)) {
        tabledata <- custom_list_split[[i]]
        summary_tables[[i]] <- glue::glue("<div><a href='?show=gene&query={tabledata$id}' target='_blank'><h3>{tabledata$id}</a>: {tabledata$protein_name}</h3></div>
                                          <div><b>Uniprot ID: </b><a href='https://www.uniprot.org/uniprot/{tabledata$uniprot_id}' target='_blank'>{tabledata$uniprot_id}</a></div>
                                          <div><b>Enzyme Commission: </b><a href='https://enzyme.expasy.org/EC/{tabledata$ec}' target='_blank'>{tabledata$ec}</a></div>
                                          <div><b>Protein Mass: </b>{tabledata$mass} kDa</div>
                                          <div><b>Description: </b></div>
                                          <div><p>{tabledata$function_cc}</p></div>
                                          ")
      }
      return(dplyr::bind_rows(summary_tables) %>%
               tidyr::unite("text", dplyr::everything(), sep = " ")
      )
    }

    valid_summaries <-
      custom_list_split %>%
      tab_fun_gene() %>%
      dplyr::pull(text) %>%
      htmltools::HTML()
  }

  return(valid_summaries)
}

#' Make Protein Sequence
#'
#' The make_protein_sequence function takes a gene as an input and returns its protein sequence.
#'
#' @param input A list containing a content variable.
#' @param var Variable that determines which text is returned
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_protein_sequence(input = list(type = "gene", content = c("ROCK1")))
#' make_protein_sequence(input = list(type = "gene", content = c("ROCK1", "ROCK2")))
make_protein_sequence <- function(input = list(),
                                  ...) {
  if (is.null(input$content)) {
    return (NULL)
  }

  valid_sequence <-
    get_data_object(object_name = input$content,
                    dataset_name = "universal_proteins",
                    pivotwider = TRUE) %>%
    dplyr::pull(sequence)

  return(valid_sequence)
}

#' Make Cell Summary
#'
#' The make_summary_cell function takes a cell as an input and returns summary text about it.
#'
#' @param input A list containing a content variable.
#' @param var Variable that determines which text is returned
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_summary_cell(input = list(content = "HEPG2"), var = "cell_line")
#' make_summary_cell(input = list(content = "HEPG2"), var = "lineage_subtype")
make_summary_cell <- function(input = list(),
                              var = "cell_line") {
  # if (is.null(input$content)) {
  #   return (NULL)
  # }
  # cell_summary_var <-
  #   data_cell_expression_names %>%
  #   dplyr::filter(cell_line %in% input$content) %>%
  #   dplyr::pull(var) #any column name
  # return(cell_summary_var)
}

#' Make Lineage Summary
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_summary_lineage(input = list(content = "Cervix"))
make_summary_lineage <- function(input = list(),
                                 var = "cell_line") { #default so no error if empty, but this pulls the var out of the df
  # if (is.null(input$query)) {
  #   return (NULL)
  # }
  # cell_lineage_var <-
  #   data_cell_expression_names %>%
  #   dplyr::filter(lineage == input$query) %>%
  #   dplyr::pull(var) #any column name
  # return(cell_lineage_var)
}

#' Make cell_osaurus Summary
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_summary_cellosaurus(input = list(content = "HEPG2"), var = "SX") #sex
#' make_summary_cellosaurus(input = list(content = "HEPG2"), var = "AG") #age
#' make_summary_cellosaurus(input = list(content = "HEPG2"), var = "CC") %>% lit_linkr(data_universal_gene_summary = universal_gene_summary)
#' make_summary_cellosaurus(input = list(content = "HEPG2"), var = "ATCC_url") #url
make_summary_cellosaurus <- function(input = list(),
                                     var = "ID") {
  # cell_var <-
  #   data_cell_osaurus %>%
  #   dplyr::filter(name %in% input$content) %>%
  #   dplyr::select(var) %>%
  #   dplyr::pull()
  #
  # return(cell_var)
}

#' Get Essential
#'
#' Gets number of cells lines a gene is essential within
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' get_essential(input = list(type = "gene", content = "ROCK1"))
#' get_essential(input = list(type = "gene", content = "ROCK2"))
#' get_essential(input = list(type = "gene", content = c("ROCK2", "ROCK1")))
#' get_essential(input = list(type = "gene", subtype = "custom_gene_list", content = c("RDX", "ROCK2", "DTX3L", "MSN", "SORL1", "EZR")))
#' get_essential(input = list(type = "cell", content = "HEPG2"))
#' get_essential(input = list(type = "cell", content = c("HEPG2", "HUH7")))
#' get_essential(input = list(type = "compound", content = "aspirin"))
get_essential <- function(input = list(),
                          ...) {
  if(input$type == "gene") {
    #new pre-computed
    essential <-
      get_data_object(object_name = input$content,
                      dataset_name = "gene_essential",
                      pivotwider = TRUE)
    if(length(input$content) >1){
      #instead of reconstructing it, you could simply grab them, and str_c(). what would be the best (most concise) format?
      #what about counting unique? (data cannot do that as presented)...
      sentence <-
        essential %>%
        dplyr::arrange(desc(n)) %>%
        glue::glue_data("{id} is essential in {n} {cell_line}") %>%
        stringr::str_c(collapse = ", ")
    } else {
      #length == 1
      sentence <-
        essential %>%
        dplyr::pull(sentence)
    }
  # } else if(input$type == "compound") {
    #ibid
    #when you start this up again, the data are ready to be generated in ddh-data, and then added to a query object (not yet generated), and then queried here using same logic
  # } else if(input$type == "cell") {
    #ibid
  # } else {
  #   stop("uh oh")
  }
  return(sentence)
}

#' Make Compound Summary
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_summary_compound(input = list(content = "aspirin"), var = "name")
make_summary_compound <- function(input = list(),
                                  var = "name") { #default so no error if empty, but this pulls the var out of the df
  # if (is.null(input$query)) {
  #   return (NULL)
  # }
  # compound_summary_var <-
  #   data_universal_prism_meta %>%
  #   dplyr::filter(name == input$content) %>%
  #   dplyr::pull(var)#any column name
  # return(compound_summary_var)
}

#' Make Metabolite Summary
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_summary_metabolite(input = list(content = "citrate"), var = "name")
make_summary_metabolite <- function(input = list(),
                                    var = "name") { #default so no error if empty, but this pulls the var out of the df
  # if (is.null(input$query)) {
  #   return (NULL)
  # }
  # compound_summary_var <-
  #   data_compound_hmdb_names %>%
  #   dplyr::filter(name == input$content) %>% #are names good enough at matching?
  #   dplyr::pull(var) #any column name
  # return(compound_summary_var)
}

#' Make MOA Summary
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_summary_moa(input = list(content = "cyclooxygenase inhibitor", content = "cyclooxygenase inhibitor"), var = "name")
#' make_summary_moa(input = list(content = "cyclooxygenase inhibitor", content = "cyclooxygenase inhibitor"), var = "moa")
make_summary_moa <- function(input = list(),
                             var = "name",
                             ...) {
  # if (is.null(input$query)) {
  #   return (NULL)
  # }
  # moa_summary_var <-
  #   data_universal_prism_meta %>%
  #   dplyr::filter(moa == input$content) %>%
  #   dplyr::pull(var) %>%
  #   stringr::str_c(collapse = ", ")
  #
  # return(moa_summary_var)
}

#' Make Compound List Summary
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_summary_compound_list(input = list(content = c("aspirin", "carprofen")))
make_summary_compound_list <- function(input = list(),
                                       ...) {
  # if (is.null(input$content)) {
  #   return (NULL)
  # }
  # # Filter out invalid symbols for when a user edits "custom_gene_list" content parameter
  # valid_compound_name <-
  #   data_universal_prism_meta %>%
  #   dplyr::filter(name %in% input$content) %>%
  #   dplyr::pull(name) %>%
  #   stringr::str_c(collapse = ", ")
  # return(valid_compound_name)
}

