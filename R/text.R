# Gene Summary
summary_gene <- function(summary_table = gene_summary,
                         input = list(),
                         var = "approved_symbol") { #default so no error if empty, but this pulls the var out of the df
  if (is.null(input$content)) {
    return (NULL)
  }
  gene_summary_var <- summary_table %>%
    dplyr::filter(approved_symbol == input$content) %>%
    dplyr::pull(var) #any column name
  return(gene_summary_var)
}
# summary_gene(input = list(content = "ROCK1"), var = "approved_symbol")
# summary_gene(input = list(content = "ROCK1"), var = "aka")
# summary_gene(summary_table = gene_location, input = list(content = "ROCK1"), var = "cds_length")

# Pathway Summary
summary_pathway <- function(summary_table = pathways, input = list(), var = "pathway") {
  if (is.null(input$query)) {
    return (NULL)
  }
  if (var == "data") {
    pathway_summary_var <- summary_table %>%
      dplyr::filter(go == input$query) %>%
      unnest(data) %>%
      dplyr::pull(gene) %>%
      str_c(collapse = ", ")
  } else {
    pathway_summary_var <- summary_table %>%
      dplyr::filter(go == input$query) %>%
      dplyr::pull(var)
  }
  return(pathway_summary_var)
}
#summary_pathway(input = list(query = "1902965"), var = "data")

# Gene list Summary
summary_list <- function(summary_gene = gene_summary,
                         location_gene = gene_location,
                         summary_cell = expression_meta,
                         cell_meta = cellosaurus,
                         summary_len = 40, # number of WORDS
                         input = list()) {
  if (is.null(input$content)) {
    return (NULL)
  }

  if(input$type == "gene") {
    custom_list <- summary_gene %>%
      dplyr::filter(approved_symbol %in% input$content) %>%
      left_join(location_gene, by = "approved_symbol")

    custom_list[custom_list == "NA"] <- NA
    custom_list[custom_list == ""] <- NA
    custom_list[is.na(custom_list)] <- "No info."

    if(length(input$content) == 1) {

      valid_summaries <- glue::glue("<div><h3>{custom_list$approved_symbol}: {custom_list$approved_name}</h3></div>
                                    <div><b>Entrez ID: </b><a href='https://www.ncbi.nlm.nih.gov/gene/?term={custom_list$ncbi_gene_id}' target='_blank'>{custom_list$ncbi_gene_id}</a></div>
                                    <div><b>ENSEMBL ID: </b>{custom_list$ensembl_gene_id}</div>
                                    <div><b>Chromosome: </b>{custom_list$chromosome}</div>
                                    <div><b>Coding Sequence Length: </b>{custom_list$cds_length} bp</div>
                                    <div><b>Aka: </b>{custom_list$aka}</div>
                                    <div><b>Description</b></div>
                                    <div><p>{custom_list$entrez_summary}</p></div>
                                    ") %>%
        HTML()
    }
    else {
      custom_list <- custom_list %>%
        mutate(entrez_summary = ifelse(str_count(entrez_summary) >= summary_len,
                                       paste0(gsub(paste0("^((\\w+\\W+){", summary_len, "}\\w+).*$"), "\\1", entrez_summary), " ..."),
                                       entrez_summary)
        )

      custom_list_split <- split(custom_list, custom_list$approved_symbol)

      tab_fun_gene <- function(custom_list_split){
        summary_tables <- list()
        for (i in names(custom_list_split)){
          tabledata <- custom_list_split[[i]]
          summary_tables[[i]] <- glue::glue("<div><a href='?show=gene&query={tabledata$approved_symbol}' target='_blank'><h3>{tabledata$approved_symbol}</a>: {tabledata$approved_name}</h3></div>
                                            <div><b>Entrez ID: </b><a href='https://www.ncbi.nlm.nih.gov/gene/?term={tabledata$ncbi_gene_id}' target='_blank'>{tabledata$ncbi_gene_id}</a></div>
                                            <div><b>ENSEMBL ID: </b>{tabledata$ensembl_gene_id}</div>
                                            <div><b>Description</b></div>
                                            <div><p>{tabledata$entrez_summary}</p></div>
                                            ")
        }
        return(bind_rows(summary_tables) %>%
                 tidyr::unite("text", everything(), sep = " ")
        )
      }

      valid_summaries <- custom_list_split %>%
        tab_fun_gene() %>%
        pull(text) %>%
        HTML()
    }

  } else if (input$type == "cell") {
    custom_list <- summary_cell %>%
      dplyr::filter(cell_line %in% input$content) %>%
      dplyr::select(cell_line, lineage, lineage_subtype, age, sex) %>%
      left_join(cell_meta %>%
                  dplyr::select(name, CC) %>%
                  dplyr::rename(cell_line = name),
                by = c("cell_line")) %>%
      replace_na(list(CC = "NA"))

    custom_list[custom_list == "NA"] <- NA
    custom_list[custom_list == ""] <- NA
    custom_list[is.na(custom_list)] <- "No info."

    if(length(input$content) == 1) {

      valid_summaries <- glue::glue("<div><h3>{custom_list$cell_line} ({custom_list$lineage})</h3></div>
                                    <div><b>Lineage: </b>{custom_list$lineage}</div>
                                    <div><b>Lineage subtype: </b>{custom_list$lineage_subtype}</div>
                                    <div><b>Age: </b>{custom_list$age}</div>
                                    <div><b>Sex: </b>{custom_list$sex}</div>
                                    <div><b>Description</b></div>
                                    <div><p>{custom_list$CC}</p></div>
                                    ") %>%
        HTML()
    }
    else {
      custom_list <- custom_list %>%
        mutate(CC = ifelse(str_count(CC) >= summary_len,
                           paste0(gsub(paste0("^((\\w+\\W+){", summary_len, "}\\w+).*$"), "\\1", CC), " ..."),
                           CC)
        )

      custom_list_split <- split(custom_list, custom_list$cell_line)

      tab_fun_cell <- function(custom_list_split){
        summary_tables <- list()
        for (i in names(custom_list_split)){
          tabledata <- custom_list_split[[i]]
          summary_tables[[i]] <- glue::glue("<div><a href='?show=cell&query={tabledata$cell_line}' target='_blank'><h3>{tabledata$cell_line} ({tabledata$lineage})</a></h3></div>
                                            <div><b>Lineage: </b>{tabledata$lineage}</div>
                                            <div><b>Lineage subtype: </b>{tabledata$lineage_subtype}</div>
                                            <div><b>Description</b></div>
                                            <div><p>{tabledata$CC}</p></div>
                                            ")
        }
        return(bind_rows(summary_tables) %>%
                 tidyr::unite("text", everything(), sep = " ")
        )
      }

      valid_summaries <- custom_list_split %>%
        tab_fun_cell() %>%
        pull(text) %>%
        HTML()
    }
  }

  return(valid_summaries)
}

#summary_list(input = list(type = "gene", content = c("ROCK1", "ROCK2")))
#summary_list(input = list(type = "cell", content = c("HEL", "HEPG2")))

# protein summary
#this is a master function to pull data out of the proteins df
summary_protein <- function(summary_table = proteins, input = list(), var = "gene_name") { #default so no error if empty, but this pulls the var out of the df
  if (is.null(input$query)) {
    return (NULL)
  }
  protein_summary_var <- summary_table %>%
    dplyr::filter(gene_name %in% input$query) %>%
    dplyr::pull(var) #any column name
  return(protein_summary_var)
}
#summary_protein(input = list(query = "ROCK1"), var = "gene_name")
#summary_protein(input = list(query = "ROCK2"), var = "protein_name")
#summary_protein(summary_table = proteins, input = list(query = "ROCK1"), var = "sequence")

# cell summary
summary_cell <- function(summary_table = expression_names,
                         input = list(),
                         var = "cell_line") { #default so no error if empty, but this pulls the var out of the df
  if (is.null(input$content)) {
    return (NULL)
  }
  cell_summary_var <-
    summary_table %>%
    dplyr::filter(cell_line %in% input$content) %>%
    dplyr::pull(var) #any column name
  return(cell_summary_var)
}
#summary_cell(input = list(cell_line = "HEPG2"), var = "cell_line")
#summary_cell(input = list(cell_line = "HEPG2"), var = "lineage_subtype")

# cell lineage summary
summary_lineage <- function(summary_table = expression_names,
                            input = list(),
                            var = "cell_line") { #default so no error if empty, but this pulls the var out of the df
  if (is.null(input$query)) {
    return (NULL)
  }
  cell_lineage_var <-
    summary_table %>%
    dplyr::filter(lineage == input$query) %>%
    dplyr::pull(var) #any column name
  return(cell_lineage_var)
}
#summary_lineage(input = list(query = "Cervix"))

### cellosaurus
summary_cellosaurus <- function(summary_table = cellosaurus,
                                key_table = cellosaurus_key,
                                input = list(),
                                var = "ID") {
  cell_var <-
    summary_table %>%
    dplyr::filter(name %in% input$content) %>%
    dplyr::select(var) %>%
    pull()

  return(cell_var)
}
#summary_cellosaurus(input = list(content = "HEPG2"), var = "SX") #sex
#summary_cellosaurus(input = list(content = "HEPG2"), var = "AG") #age
#summary_cellosaurus(input = list(content = "HEPG2"), var = "CC") %>% lit_linkr(summary_table = gene_summary) #
#summary_cellosaurus(input = list(content = "HEPG2"), var = "ATCC_url") #url

#gets number of cells lines a gene is essential within
get_essential <- function(achilles_data = achilles_long,
                          prism_data = prism_long,
                          input = list()) {
  if(input$type == "gene") {
    essential <-
      achilles_data %>%
      dplyr::filter(gene %in% input$content) %>%
      dplyr::group_by(gene) %>%
      dplyr::summarize(n = sum(dep_score < -1, na.rm = TRUE)) %>%
      dplyr::mutate(cell_line = case_when(
        n == 1 ~ "cell line",
        TRUE ~ "cell lines"
      )) %>%
      dplyr::arrange(desc(n)) %>%
      glue::glue_data("{gene} is essential in {n} {cell_line}") %>%
      stringr::str_c(collapse = ", ")
  } else if(input$type == "compound") {
    essential <-
      prism_data %>%
      dplyr::filter(name %in% input$content) %>%
      dplyr::group_by(name) %>%
      dplyr::summarize(n = sum(log2fc < -1, na.rm = TRUE)) %>%
      dplyr::mutate(cell_line = case_when(
        n == 1 ~ "cell line",
        TRUE ~ "cell lines"
      )) %>%
      dplyr::arrange(desc(n)) %>%
      glue::glue_data("{name} is toxic in {n} {cell_line}") %>%
      stringr::str_c(collapse = ", ")
  } else if(input$type == "cell") {
    essential <-
      achilles_data %>%
      dplyr::left_join(expression_names, by = "X1") %>%
      dplyr::select(cell_line, gene, dep_score) %>%
      dplyr::filter(cell_line %in% input$content) %>%
      dplyr::group_by(cell_line) %>%
      dplyr::summarize(n = sum(dep_score < -1, na.rm = TRUE)) %>%
      dplyr::mutate(gene = case_when(
        n == 1 ~ "gene",
        TRUE ~ "genes"
      )) %>%
      dplyr::arrange(desc(n)) %>%
      glue::glue_data("{cell_line} has {n} essential {gene}") %>%
      stringr::str_c(collapse = ", ")
  } else {
    stop("uh oh")
  }
  return(essential)
}

#get_essential(input = list(type = "gene", content = "ROCK1"))
#get_essential(input = list(type = "gene", content = "ROCK2"))
#get_essential(input = list(type = "gene", content = c("ROCK2", "ROCK1")))
#get_essential(input = list(type = "gene", query = "custom_gene_list", content = c("RDX", "ROCK2", "DTX3L", "MSN", "SORL1", "EZR")))
#get_essential(input = list(type = "cell", content = "HEPG2"))
#get_essential(input = list(type = "cell", content = c("HEPG2", "HUH7")))
#get_essential(input = list(type = "compound", content = "aspirin"))

#this group contains master functions for returning text from summary tables
# Compound Summary
summary_compound <- function(summary_table = prism_meta,
                             input = list(),
                             var = "name") { #default so no error if empty, but this pulls the var out of the df
  if (is.null(input$query)) {
    return (NULL)
  }
  compound_summary_var <-
    summary_table %>%
    dplyr::filter(name == input$content) %>%
    dplyr::pull(var)#any column name
  return(compound_summary_var)
}
# summary_compound(input = list(query = "aspirin"), var = "name")

summary_metabolite <- function(summary_table = hmdb_names,
                               input = list(),
                               var = "name") { #default so no error if empty, but this pulls the var out of the df
  if (is.null(input$query)) {
    return (NULL)
  }
  compound_summary_var <-
    summary_table %>%
    dplyr::filter(cid == input$content) %>%
    dplyr::pull(var)#any column name
  return(compound_summary_var)
}

# MOA Summary
summary_moa <- function(summary_table = prism_meta,
                        input = list(),
                        var = "name") {
  if (is.null(input$query)) {
    return (NULL)
  }
  moa_summary_var <-
    summary_table %>%
    dplyr::filter(moa == input$content) %>%
    dplyr::pull(var) %>%
    str_c(collapse = ", ")

  return(moa_summary_var)
}
#summary_moa(input = list(query = "cyclooxygenase inhibitor"), var = "name")
#summary_moa(input = list(query = "cyclooxygenase inhibitor"), var = "moa")

# Compound list Summary
summary_compound_list <- function(summary_table = prism_meta,
                                  input = list()) {
  if (is.null(input$content)) {
    return (NULL)
  }
  # Filter out invalid symbols for when a user edits "custom_gene_list" content parameter
  valid_compound_name <-
    summary_table %>%
    dplyr::filter(name %in% input$content) %>%
    dplyr::pull(name) %>%
    str_c(collapse = ", ")
  return(valid_compound_name)
}
# summary_compound_list(input = list(content = c("aspirin", "carprofen")))

