# SETUP GRAPH ----------------------------------------------------------------------
#' Setup graph parameters
#'
#' The overall purpose of this function is to create an object that can be used to generate a network graph by filtering the query object. A list is returned that has four elements.
#'
#' @importFrom magrittr %>%
#'
#' @export
setup_graph <- function(setup_input = list(), #changed name here to prevent var naming overlap for nested funs()
                        setup_threshold,
                        setup_corr_type,
                        setup_card) {

  pathway_ids <- ddh::get_content("universal_pathways", dataset = TRUE) %>%
    dplyr::distinct(gs_id, gs_name) %>%
    dplyr::mutate_all(as.character) %>%
    dplyr::mutate(set = gsub("_.*", "", gs_name)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(gs_name = gsub(paste0(set, "_"), "", gs_name)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(gs_name = gsub("_", " ", gs_name))

  if (setup_input$query %in% pathway_ids$gs_id) { # pathways
    # get master data object, which has all pathway + related
    dep_network_master <-
      get_data_object(object_name = setup_input$query,
                      dataset_name = "setup_graph") %>%
      tidyr::separate(col = "key", into = c("type", "rank"), sep = "_") %>%
      dplyr::mutate(rank = as.numeric(rank)) %>%
      dplyr::distinct(id, data_set, rank, value) # remove redundancy

    threshold_pathways <-
      dep_network_master %>%
      dplyr::filter(id %in% setup_input$query,
                    rank <= setup_threshold) %>%
      dplyr::pull(., value) %>%
      unique()

    threshold_pathways <- c(setup_input$query, threshold_pathways)

    #this makes dep_network_table tibble to generate the network graph
    dep_network_table <-
      dep_network_master %>%
      dplyr::filter(id %in% threshold_pathways,
                    rank <= setup_threshold)

    setup_object <- list(dep_network_table = dep_network_table,
                         threshold_pathways = threshold_pathways, # pathways used to create graph
                         pathway_ids = pathway_ids)

  } else { # genes
    #set corr_filter from corr_type
    if(setup_corr_type == "both"){corr_filter = c("positive", "negative")} else {corr_filter = setup_corr_type}

    #filter for cards
    if (setup_card == TRUE & length(setup_input$content) > 5) {
      setup_input$content <- sample(setup_input$content, 5)
    }

    #get master data object, which has all genes + related
    #some redundancy b/c if i'm in your top 10, and you're in mine, I fetch you twice when two feather objects come together
    dep_network_master <-
      get_data_object(object_name = setup_input$content,
                      dataset_name = "setup_graph") %>%
      dplyr::filter(stringr::str_detect(.$key, "negative") | stringr::str_detect(.$key, "positive")) %>% #need this so approved_name doesn't get separated in next step
      tidyr::separate(col = "key", into = c("type", "rank"), sep = "_") %>%
      dplyr::mutate(rank = as.numeric(rank)) %>%
      dplyr::distinct(id, data_set, type, rank, value) #remove redundancy

    #this is the master threshold gene vec to use for filtering and factor grouping
    threshold_genes_pos <- NULL
    threshold_genes_neg <- NULL

    #set positive and/or negative thresholds
    if("positive" %in% corr_filter){
      threshold_genes_pos <-
        dep_network_master %>%
        dplyr::filter(id %in% setup_input$content,
                      type == "positive",
                      rank <= setup_threshold) %>%
        dplyr::pull(., value) %>%
        unique()
    }
    #need to also do for negative, in case we need them
    if("negative" %in% corr_filter){
      threshold_genes_neg <-
        dep_network_master %>%
        dplyr::filter(id %in% setup_input$content,
                      type == "negative",
                      rank <= setup_threshold) %>%
        dplyr::pull(., value) %>%
        unique()
    }

    #this next step is key: either pull top n genes for single query, or pull query genes from multi-gene query
    if(length(setup_input$content) == 1){
      #get single gene threshold vec, so I can use this to filter my id col for the dep_network_table
      threshold_genes <- c(setup_input$content, threshold_genes_pos, threshold_genes_neg)
    } else {
      #this only keeps input$content in id, and drops related
      threshold_genes <- setup_input$content
    }

    #this makes dep_network_table tibble to generate the network graph
    dep_network_table <-
      dep_network_master %>%
      dplyr::filter(id %in% threshold_genes,
                    type %in% corr_filter,
                    rank <= setup_threshold)

    setup_object <- list(dep_network_table = dep_network_table, #full dataset
                         threshold_genes_pos = threshold_genes_pos, #req'd for graph factors/labels
                         threshold_genes_neg = threshold_genes_neg, #req'd for graph factors/labels
                         threshold_genes = threshold_genes) #genes used to create graph
  }
  return(setup_object)
}

## MAKE GRAPH ----------------------------------------------------------------------
#' Create network graph visualization using visNetwork
#'
#' This function takes in dependency correlations and a gene query list to then output a dependency network graph
#' visualization containing the top/bottom threshold for each of the top/bottom threshold of the gene query list
#' using visNetwork.
#'
#' @param input A list containing character vector of gene_symbols used to create network graph
#' @param threshold A numerical representing the number of genes to pull from top and bottom tables
#' @param deg A numerical representing the minimum number of connections for a gene to be connected to the network
#' @param corr_type A string that describes what type of correlations to include, options are: "positive", "negative", or "both"
#' @param displayHeight Default to "90vh". The height of the network in pixels ("500px"), as a percentage (100 percent), or as a percentage of the viewport ("70vh", where 70 represents 70 percent of the viewport)
#' @param displayWidth Default to 100 percent. The width of the network in pixels ("500px"), as a percentage (100 percent), or as a percentage of the viewport ("70vh", where 70 represents 70 percent of the viewport)
#' @param tooltipLink Boolean to denote whether or not to include a link in the tooltip for a gene. Default to false.
#'
#' @return Outputs a complete network graph. If an error is thrown, then will return an empty graph.
#'
#' @importFrom magrittr %>%
#' @import visNetwork
#'
#' @export
#' @examples
#' make_graph(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_graph(input = list(type = "gene", query = 'ROCK1', content = "ROCK1"), card = TRUE)
#' make_graph(input = list(type = "gene", query = 'ROCK1', content = "ROCK1"), corr_type = "negative")
#' make_graph(input = list(type = "gene", query = 'ROCK1', content = "ROCK1"), corr_type = "both")
#' make_graph(input = list(type = "gene", query = 'ROCK1', content = c("ROCK1", "ROCK2")))
#' make_graph(input = list(type = "pathway", query = "16769"))
#' make_graph(input = list(type = "pathway", query = "16769"), tooltipLink = TRUE)
#' make_graph(input = list(type = "gene", query = 'DTX3L', content = "DTX3L"), corr_type = "negative") # disconnected query gene
#' make_graph(input = list(type = "compound", query = 'aspirin', content = "aspirin"), corr_type = "negative")
make_graph <- function(input = list(),
                       threshold = 10,
                       deg = 2,
                       cell_line_var = "dependency",
                       # data_cell_line_dep_sim = cell_dependency_sim,
                       # data_cell_line_exp_sim = cell_expression_sim,
                       bonferroni_cutoff = 0.05,
                       corr_type = "positive",
                       displayHeight = '90vh',
                       displayWidth = '100%',
                       tooltipLink = FALSE,
                       card = FALSE) {
  make_graph_raw <- function() {
    #set color schemes
    ddh::load_ddh_colors()
    if (input$type == "gene" | input$type == "pathway") {
      queryColor <- color_set_gene_alpha[2]
    }  else if (input$type == "cell"){
      queryColor <- color_set_cell_alpha[2]
    } else if(input$type == "compound") {
      queryColor <- color_set_compound_alpha[2]
    } else {
      stop("declare your type")
    }

    #build graph using setup
    setup_graph_list <- setup_graph(setup_input = input, #changed name here to prevent var naming overlap for nested funs()
                                    setup_threshold = threshold,
                                    setup_corr_type = corr_type,
                                    setup_card = card) #adds the filter logic to the data

    disconnected <- FALSE

    if (input$type == "pathway") {
      corr_type <- "pathway"
      # make graph
      graph_network <-
        setup_graph_list$dep_network_table %>%
        dplyr::rename(from = id, to = value) %>%
        tidygraph::as_tbl_graph()

      group_var <- c("Query", "Co-essential", "Connected")

      # add data to nodes
      # query_name <- setup_graph_list$pathway_ids[which(setup_graph_list$pathway_ids$gs_id == input$query), 2]$gs_name
      # threshold_pathways_name <- setup_graph_list$pathway_ids[which(setup_graph_list$pathway_ids$gs_id %in% setup_graph_list$threshold_pathways), 2]$gs_name

      nodes <- # active by default
        graph_network %>%
        tidygraph::as_tibble() %>%
        tibble::rowid_to_column("id") %>%
        dplyr::mutate(degree = igraph::degree(graph_network),
                      group = dplyr::case_when(name %in% input$query ~ "Query",
                                               name %in% setup_graph_list$threshold_pathways ~ "Co-essential",
                                               TRUE ~ "Connected"),
                      group = forcats::as_factor(group),
                      group = forcats::fct_relevel(group, group_var)) %>%
        dplyr::arrange(group)

      links <-
        graph_network %>%
        tidygraph::activate(edges) %>%
        tidygraph::as_tibble()

      # determine the nodes that have at least the minimum degree
      nodes_filtered <-
        nodes %>%
        dplyr::filter(degree >= deg) %>%  #input$degree
        as.data.frame()

      # filter the edge list to contain only links to or from the nodes that have the minimum or more degree
      links_filtered <-
        links %>%
        dplyr::filter(to %in% nodes_filtered$id & from %in% nodes_filtered$id) %>%
        as.data.frame()

      # re-adjust the from and to values to reflect the new positions of nodes in the filtered nodes list
      links_filtered$from <- match(links_filtered$from, nodes_filtered$id) - 1
      links_filtered$to <- match(links_filtered$to, nodes_filtered$id) - 1

      #check to see if setting degree removed all links; if so, then throws error, so this fills a dummy links_filtered df to plot only nodes
      if(nrow(links_filtered) == 0) {
        links_filtered <- dplyr::tibble("from" = -1, "to" = -1, "cc" = 1, "origin" = "coessential")
      }

      # shift id values properly to work with visNetwork
      nodes_filtered <-
        nodes_filtered %>%
        dplyr::mutate(id=0:(nrow(nodes_filtered)-1))

      # recreate the network using filtered edges to use degree function
      graph_network_filtered <-
        igraph::graph_from_data_frame(links_filtered, directed = F) %>%
        igraph::simplify()
      degVec <- igraph::degree(graph_network_filtered)

      # make node size a function of their degree in the current network
      if(names(degVec)[1] == "-1"){ # if there are no links remaining then assign degree of each node to 0
        nodes_filtered <-
          nodes_filtered %>%
          dplyr::mutate(value = 0)
      } else if(length(input$content > 1) & length(degVec) != dim(nodes_filtered)[1]){ # handle cases where some query genes are disconnected
        genesWithConnections <-
          degVec %>%
          names() %>%
          as.numeric()
        nodes_filtered <-
          nodes_filtered %>%
          tibble::add_column(value = 0)
        for(gene in 1:dim(nodes_filtered)[1]){
          if(nodes_filtered[gene,"id"] %in% genesWithConnections){
            nodes_filtered[gene,"value"] <- degVec[nodes_filtered[gene,"id"] %>% toString()]
          }
        }
      } else {
        nodes_filtered <-
          nodes_filtered[degVec %>% names() %>% as.numeric() +1,] %>%
          dplyr::mutate(value = degVec)
      }

      # make sure query gene is at the start so legend shows proper order
      nodes_filtered <-
        nodes_filtered %>%
        dplyr::arrange(id)

      }
    else { # genes
      #make graph
      graph_network <-
        setup_graph_list$dep_network_table %>%
        dplyr::rename(from = id, to = value) %>%
        tidygraph::as_tbl_graph()

      #make groups for fct_relevel below
      corr_var <- switch(corr_type,
                         positive = "Positive",
                         negative = "Negative",
                         both = c("Positive", "Negative")
      )

      if(length(input$content) == 1){connected_var <- "Connected"} else {connected_var <- NULL}
      group_var <- c("Query", corr_var, connected_var)

      #add data to nodes
      nodes <- #active by default
        graph_network %>%
        tidygraph::as_tibble() %>%
        tibble::rowid_to_column("id") %>%
        dplyr::mutate(degree = igraph::degree(graph_network),
                      #name is the list of nodes
                      group = dplyr::case_when(name %in% input$content == TRUE ~ "Query",
                                               name %in% setup_graph_list$threshold_genes_pos == TRUE ~ "Positive",
                                               name %in% setup_graph_list$threshold_genes_neg == TRUE ~ "Negative",
                                               TRUE ~ "Connected"),
                      group = forcats::as_factor(group),
                      group = forcats::fct_relevel(group, group_var)) %>%
        dplyr::arrange(group)

      links <-
        graph_network %>%
        tidygraph::activate(edges) %>% # %E>%
        tidygraph::as_tibble()

      # determine the nodes that have at least the minimum degree
      nodes_filtered <-
        nodes %>%
        dplyr::filter(degree >= deg) %>%  #input$degree
        as.data.frame

      # filter the edge list to contain only links to or from the nodes that have the minimum or more degree
      links_filtered <-
        links %>%
        dplyr::filter(to %in% nodes_filtered$id & from %in% nodes_filtered$id) %>%
        as.data.frame

      # re-adjust the from and to values to reflect the new positions of nodes in the filtered nodes list
      links_filtered$from <- match(links_filtered$from, nodes_filtered$id) - 1
      links_filtered$to <- match(links_filtered$to, nodes_filtered$id) - 1

      #check to see if setting degree removed all links; if so, then throws error, so this fills a dummy links_filtered df to plot only nodes
      if(nrow(links_filtered) == 0) {
        if(corr_type == "Negative"){
          links_filtered <- dplyr::tibble("from" = -1, "to" = -1, "r2" = 1, "origin" = "neg")
        }    else{
          links_filtered <- dplyr::tibble("from" = -1, "to" = -1, "r2" = 1, "origin" = "pos")}
      }

      # shift id values properly to work with visNetwork
      nodes_filtered <-
        nodes_filtered %>%
        dplyr::mutate(id=0:(nrow(nodes_filtered)-1))

      # recreate the network using filtered edges to use degree function
      graph_network_filtered <-
        igraph::graph_from_data_frame(links_filtered, directed = F) %>%
        igraph::simplify()
      degVec <- igraph::degree(graph_network_filtered)

      # make node size a function of their degree in the current network
      if(names(degVec)[1] == "-1"){ # if there are no links remaining then assign degree of each node to 0
        nodes_filtered <-
          nodes_filtered %>%
          dplyr::mutate(value = 0)
      } else if(length(input$content > 1) & length(degVec) != dim(nodes_filtered)[1]){ # handle cases where some query genes are disconnected
        genesWithConnections <-
          degVec %>%
          names() %>%
          as.numeric()
        nodes_filtered <-
          nodes_filtered %>%
          tibble::add_column(value = 0)
        for(gene in 1:dim(nodes_filtered)[1]){
          if(nodes_filtered[gene,"id"] %in% genesWithConnections){
            nodes_filtered[gene,"value"] <- degVec[nodes_filtered[gene,"id"] %>% toString()]
          }
        }
      } else {
        nodes_filtered <-
          nodes_filtered[degVec %>% names() %>% as.numeric() +1,] %>%
          dplyr::mutate(value = degVec)
      }

      # make sure query gene is at the start so legend shows proper order
      nodes_filtered <-
        nodes_filtered %>%
        dplyr::arrange(id)

      #check to see if query gene is missing; if so, then adds a dummy so it shows up on graph, but disconnected
      # disconnected <- FALSE
      #could simplify this if/else with common input name, b/c only here to get input$content vs. input$content
      # if(sum(str_detect(nodes_filtered$group, "Query")) == 0){
      #   dummy <- tibble("id" = max(nodes_filtered$id) + 1, "name" = input$content, "degree" = 1, "group" = "Query", "value" = 0)
      #   nodes_filtered <- bind_rows(dummy, nodes_filtered)
      #   disconnected <- T
      # } else {
      #   stop("declare your type")
      # }

      # if(input$type == "cell") {
      #   for(cell in nodes_filtered$name){
      #     newVal <-
      #       cell_top_data %>%
      #       dplyr::filter(cell1_name == cell | cell2_name == cell)
      #
      #     # Swap cols (based on query)
      #     for(i in 1:nrow(newVal)) {
      #       if(newVal$cell2_name[i] %in% cell & !(newVal$cell1_name[i] %in% cell)) {
      #         cell1 <- newVal$cell1_name[i]
      #         cell2 <- newVal$cell2_name[i]
      #
      #         newVal$cell2_name[i] <- cell1
      #         newVal$cell1_name[i] <- cell2
      #       }
      #     }
      #
      #     newVal <-
      #       newVal %>%
      #       dplyr::pull(cell2_name)
      #
      #     if(length(newVal) == 0){
      #       nameTable <- tibble::add_row(nameTable, name = "No cell line available")# handles cases where the gene is not in the gene data_gene_summary table
      #     } else{
      #       nameTable <- tibble::add_row(nameTable, name=newVal)
      #     }
      #   }
      # } else if(input$type == "compound") {
      #   for(drug in nodes_filtered$name){
      #     newVal <- prism_meta %>%
      #       dplyr::filter(name==drug) %>%
      #       dplyr::pull(moa)
      #     if(length(newVal)==0){
      #       nameTable <- tibble::add_row(nameTable, name = "No drug MOA available")# handles cases where the gene is not in the gene data_gene_summary table
      #     } else{
      #       nameTable <- tibble::add_row(nameTable, name=newVal)
      #     }
      #   }
      # } else {
      #   stop("declare your type")
      # }

      #make df for approved_name & description of each gene nodes tibble tooltips
      gene_names <-
        get_data_object(object_name = input$content,
                        dataset_name = "setup_graph") %>%
        dplyr::filter(key == "approved_name") %>%
        dplyr::select(name = id, description = value) %>%
        dplyr::distinct(name, .keep_all = TRUE)
    }

  # add title information (tooltip that appears on hover)
  if(!tooltipLink){ # Do not form a url when just making the standalone graph while testing or making reports
    if (input$type == "pathway") {
      nodes_filtered <-
        nodes_filtered %>%
        dplyr::left_join(setup_graph_list$pathway_ids, by = c("name" = "gs_id")) %>%
        dplyr::mutate(title = paste0("<center><p>", gs_name,"<br>", set, '</p>'),
                      label = name) # this is the node name on the network
    } else {
      nodes_filtered <-
        nodes_filtered %>%
        dplyr::left_join(gene_names) %>%
        dplyr::mutate(title = paste0("<center><p>", name,"<br>", description, '</p>'),
                      label = name)
    }
  } else {
    if(input$type == "gene") {
      nodes_filtered <-
        nodes_filtered %>%
        dplyr::left_join(gene_names) %>%
        dplyr::mutate(title=paste0("<center><p>", name,"<br>",description ,'<br><a target="_blank" href="?show=gene&query=',name,'">Gene Link</a></p>'),
                      label = name )
    } else if(input$type == "pathway") {
      nodes_filtered <-
        nodes_filtered %>%
        dplyr::left_join(setup_graph_list$pathway_ids, by = c("name" = "gs_id")) %>%
        dplyr::mutate(title = paste0("<center><p>", gs_name,"<br>", set, '<br><a target="_blank" href="?show=pathway&query=', name, '">Pathway Link</a></p>'),
                      label = name) # this is the node name on the network
    } else if(input$type == "cell") {
      nodes_filtered <-
        nodes_filtered %>%
        dplyr::left_join(gene_names) %>%
        dplyr::mutate(title=paste0("<center><p>", name,"<br>",name ,'<br><a target="_blank" href="?show=cell&query=',name,'">Cell Line Link</a></p>'),
                      label = name )
    } else if(input$type == "compound") {
      nodes_filtered <-
        nodes_filtered %>%
        dplyr::left_join(gene_names) %>%
        dplyr::mutate(title=paste0("<center><p>", name,"<br>",name ,'<br><a target="_blank" href="?show=compound&query=',name,'">Drug Link</a></p>'),
                      label = name )
    } else {
      stop("declare your type")
    }
  }

  # colors used within the network
  #queryColor set above based on query type
  positiveColor <-"rgba(173, 103, 125, 0.8)"
  negativeColor <- "rgba(12, 35, 50, 0.8)"
  connectedColor <- "rgba(255, 255, 255, 0.8)" #formerly purple "rgba(84, 64, 151, 0.8)"
  borderColor <- "rgba(204, 204, 204, 0.8)" #(gray80), formerly white 255, 255, 255
  edgeColor <- "rgba(84, 84, 84, 1)"

  # Physics parameters defining the visNetwork
  iter <- 150 # number of iterations to perform of stablization before display
  gravity <- 0.5
  damping <- 0.11
  timestep <- 0.25 # reducing the timestep reduces the jitteriness of the graph and can help stabilize it
  if(disconnected){ # set up function to be called at the end of stabilization, changes the zoom to handle when query gene is disconnected and flies out of the network
    stabilizationZoomFn <- "function() {this.moveTo({scale:0.35})}"
  } else{
    stabilizationZoomFn <- "function() {}"
  }

  if(card == TRUE){
    displayHeight = "200px"
    displayWidth = "200px"
  }
  # build the network visualization
  if(corr_type == "both"){
    visNetwork::visNetwork(nodes = nodes_filtered, edges = links_filtered, width = displayWidth, height = displayHeight) %>%
      visNetwork::visOptions(highlightNearest = list(enabled = T)) %>%
      visNetwork::visGroups(groupname = "Query", color = list(background = queryColor, border =borderColor, highlight = queryColor, hover = queryColor ), shape=ifelse(cell_line_var == "dependency", "dot", "diamond"), borderWidth = 2) %>%
      visNetwork::visGroups(groupname = "Positive", color = list(background = positiveColor, border = borderColor, highlight = positiveColor, hover = positiveColor), shape=ifelse(cell_line_var == "dependency", "dot", "diamond"), borderWidth = 2) %>%
      visNetwork::visGroups(groupname = "Negative", color = list(background = negativeColor, border = borderColor, highlight = negativeColor, hover = negativeColor), shape=ifelse(cell_line_var == "dependency", "dot", "diamond"), borderWidth = 2) %>%
      visNetwork::visGroups(groupname = "Connected", color = list(background = connectedColor, border = borderColor, highlight = connectedColor, hover = connectedColor), shape=ifelse(cell_line_var == "dependency", "dot", "diamond"), borderWidth = 2) %>%
      visNetwork::visEdges(color = edgeColor, smooth = F) %>%
      visNetwork::visNodes(scaling = list(min = 10, max =20)) %>%
      visNetwork::visPhysics(barnesHut = list(damping = damping, centralGravity = gravity), timestep = timestep, stabilization = list(iterations = iter)) %>%
      visNetwork::visEvents(stabilizationIterationsDone = stabilizationZoomFn)
    # visConfigure(enabled=TRUE) # use to test out new features to add from visNetwork
  }  else if(corr_type == "positive"){
    visNetwork::visNetwork(nodes = nodes_filtered, edges = links_filtered, width = displayWidth, height = displayHeight) %>%
      visNetwork::visOptions(highlightNearest = list(enabled = T)) %>%
      visNetwork::visGroups(groupname = "Query", color = list(background = queryColor, border =borderColor, highlight = queryColor, hover = queryColor ), shape=ifelse(cell_line_var == "dependency", "dot", "diamond"), borderWidth = 2) %>%
      visNetwork::visGroups(groupname = "Positive", color = list(background = positiveColor, border = borderColor, highlight = positiveColor, hover = positiveColor), shape=ifelse(cell_line_var == "dependency", "dot", "diamond"), borderWidth = 2) %>%
      visNetwork::visGroups(groupname = "Connected", color = list(background = connectedColor, border = borderColor, highlight = connectedColor, hover = connectedColor), shape=ifelse(cell_line_var == "dependency", "dot", "diamond"), borderWidth = 2) %>%
      visNetwork::visEdges(color = edgeColor, smooth = F) %>%
      visNetwork::visNodes(scaling = list(min = 10, max =20)) %>%
      visNetwork::visPhysics(barnesHut = list(damping = damping, centralGravity = gravity), timestep = timestep, stabilization = list(iterations = iter)) %>%
      visNetwork::visEvents(stabilizationIterationsDone = stabilizationZoomFn)
  }  else if(corr_type == "negative"){
    visNetwork::visNetwork(nodes = nodes_filtered, edges = links_filtered, width = displayWidth, height = displayHeight) %>%
      visNetwork::visOptions(highlightNearest = list(enabled = T)) %>%
      visNetwork::visGroups(groupname = "Query", color = list(background = queryColor, border =borderColor, highlight = queryColor, hover = queryColor ), shape=ifelse(cell_line_var == "dependency", "dot", "diamond"), borderWidth = 2) %>%
      visNetwork::visGroups(groupname = "Negative", color = list(background = negativeColor, border = borderColor, highlight = negativeColor, hover = negativeColor), shape=ifelse(cell_line_var == "dependency", "dot", "diamond"), borderWidth = 2) %>%
      visNetwork::visGroups(groupname = "Connected", color = list(background = connectedColor, border = borderColor, highlight = connectedColor, hover = connectedColor), shape=ifelse(cell_line_var == "dependency", "dot", "diamond"), borderWidth = 2) %>%
      visNetwork::visEdges(color = edgeColor, smooth = F) %>%
      visNetwork::visNodes(scaling = list(min = 10, max = 20)) %>%
      visNetwork::visPhysics(barnesHut = list(damping = damping, centralGravity = gravity), timestep = timestep, stabilization = list(iterations = iter)) %>%
      visNetwork::visEvents(stabilizationIterationsDone = stabilizationZoomFn)
  } else if (corr_type == "pathway") {
    visNetwork::visNetwork(nodes = nodes_filtered, edges = links_filtered, width = displayWidth, height = displayHeight) %>%
      visNetwork::visOptions(highlightNearest = list(enabled = T)) %>%
      visNetwork::visGroups(groupname = "Query", color = list(background = queryColor, border =borderColor, highlight = queryColor, hover = queryColor ), shape=ifelse(cell_line_var == "dependency", "dot", "diamond"), borderWidth = 2) %>%
      visNetwork::visGroups(groupname = "Co-essential", color = list(background = positiveColor, border = borderColor, highlight = positiveColor, hover = positiveColor), shape=ifelse(cell_line_var == "dependency", "dot", "diamond"), borderWidth = 2) %>%
      visNetwork::visGroups(groupname = "Connected", color = list(background = connectedColor, border = borderColor, highlight = connectedColor, hover = connectedColor), shape=ifelse(cell_line_var == "dependency", "dot", "diamond"), borderWidth = 2) %>%
      visNetwork::visEdges(color = edgeColor, smooth = F) %>%
      visNetwork::visNodes(scaling = list(min = 10, max =20)) %>%
      visNetwork::visPhysics(barnesHut = list(damping = damping, centralGravity = gravity), timestep = timestep, stabilization = list(iterations = iter)) %>%
      visNetwork::visEvents(stabilizationIterationsDone = stabilizationZoomFn)
    }
  }
  #error handling
  tryCatch(make_graph_raw(),
           error = function(e){
             message(e)
             make_empty_graph()
           })
}

#' Bipartite Graph Graph
#'
#' \code{make_bipartite_graph} returns an image of ...
#'
#' This is a graph function that takes a gene name and returns a bipartite graph graph
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a bipartite graph graph. If an error is thrown, then will return an empty graph.
#'
#' @importFrom magrittr %>%
#' @import visNetwork
#'
#' @export
#' @examples
#' make_bipartite_graph(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_bipartite_graph(input = list(type = "gene", content = "ROCK1"), collapsed = TRUE, threshold = 10, corr_type = "positive")
#' make_bipartite_graph(input = list(type = "gene", content = c("ROCK1", "ROCK2")), collapsed = TRUE, threshold = 10, corr_type = "positive", censor = c("ADP", "Adenosine triphosphate"))
#' make_bipartite_graph(input = list(type = "compound", content = "adp"))
#' \dontrun{
#' make_bipartite_graph(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_bipartite_graph <- function(input = list(),
                                 # data_master_top_table = gene_master_top_table,
                                 # data_master_bottom_table = gene_master_bottom_table,
                                 # data_hmdb_proteins = compound_hmdb_proteins,
                                 # data_hmdb_metabolites = compound_hmdb_metabolites,
                                 # data_gene_summary = universal_gene_summary,
                                 censor = character(), #removes most common metabolites
                                 collapsed = TRUE,
                                 threshold = 10,
                                 corr_type = "positive",
                                 card = FALSE) {
  make_bipartite_graph_raw <- function() {
    #set color schemes
    ddh::load_ddh_colors()
    if(input$type == "gene") {
      #build graph using setup
      setup_graph_list <- setup_graph(setup_input = input, #changed name here to prevent var naming overlap for nested funs()
                                      setup_threshold = threshold,
                                      setup_corr_type = corr_type,
                                      setup_card = card)
      #get gene_names
      genes <- setup_graph_list$threshold_genes

      #collapsed var
      if(collapsed == TRUE) {
        collapsed_var <- "metabolite_collapsed"
      } else {
        collapsed_var <- "metabolite"
      }

      #first two cols to generate graph
      hmdb_network <-
        ddh::get_data_object(object_name = input$content,
                             dataset_name = "setup_graph") %>%
        dplyr::filter(id %in% genes,
                      key == collapsed_var) %>%
        dplyr::distinct(gene_name = id,
                        metabolite_name = value) #does a select/rename here also

      #check for nrow
      if(nrow(hmdb_network) == 0) {
        return("No associated metabolites")
      }

      #get genes that have >1 connection for filtering (so we don't have a bunch of singletons???)
      connected_genes <-
        hmdb_network %>%
        dplyr::count(gene_name, sort = TRUE) %>%
        dplyr::filter(n > 1) %>%
        dplyr::pull(gene_name)

      hmdb_filtered <-
        hmdb_network %>%
        dplyr::filter(gene_name %in% connected_genes)

      #hackish way to reorder legend, which appears based on first instance of  gene in group assignment
      hmdb_query <- #get query genes
        hmdb_filtered %>%
        dplyr::filter(gene_name %in% input$content)

      hmdb_notquery <- #get all others
        hmdb_filtered %>%
        dplyr::filter(!gene_name %in% input$content)

      hmdb_filtered <- #bind rows with query first
        hmdb_query %>%
        dplyr::bind_rows(hmdb_notquery)

      #remove common metabolites; receives character vec from shiny with pre-defined choices
      if(length(censor) > 0) {
        hmdb_filtered <-
          hmdb_filtered %>%
          dplyr::filter(!metabolite_name %in% censor)
      }
    } #extra curly due to commenting
    # } else if (input$type == "compound") {
    #   #first two cols to generate graph
    #   hmdb_network <-
    #     data_hmdb_metabolites %>%
    #     dplyr::filter(fav_metabolite %in% input$content)
    #
    #   #check for nrow
    #   if(nrow(hmdb_network) == 0) {
    #     return("No associated proteins")
    #   } else {
    #     hmdb_filtered <-
    #       hmdb_network %>%
    #       tidyr::unnest(cols = c(data)) %>%
    #       dplyr::ungroup() %>%
    #       dplyr::select(gene_name, metabolite_name) %>%
    #       dplyr::arrange(gene_name)
    #   }
    # } else {
    #   stop("delcare your type")
    # }
    #make bipartite graph
    metabolic_network <-
      igraph::graph_from_data_frame(d=hmdb_filtered,
                                    directed=FALSE)

    #Now we can create a network object from this merged information. In order to keep track of what nodes are part of each mode (metabolites or proteins) we’ll add a type column to the node data that will get a TRUE value if it is one of the proteins
    igraph::V(metabolic_network)$type <- igraph::V(metabolic_network)$name %in% hmdb_filtered$gene_name #accession
    igraph::V(metabolic_network)$group <- dplyr::case_when(
      igraph::V(metabolic_network)$name %in% input$content ~ "Query",
      igraph::V(metabolic_network)$name %in% input$content ~ "Query",
      igraph::V(metabolic_network)$name %in% hmdb_filtered$gene_name ~ "Protein",
      TRUE ~ "Metabolite") #accession

    gene_names <-
      get_data_object(object_name = input$content,
                      dataset_name = "setup_graph",
      ) %>%
      dplyr::filter(key == "approved_name") %>%
      dplyr::select(name = id, description = value) %>%
      dplyr::distinct(name, .keep_all = TRUE)

    # add title information (tooltip that appears on hover)
    igraph::V(metabolic_network)$title <- igraph::V(metabolic_network)$name %>%
      purrr::map_chr(~ glue::glue('<center>
                                <p>
                                {.x}
                                <br>
                                {if(length(gene_names %>% dplyr::filter(name == .x) %>% dplyr::pull(description)) != 0) {
                                gene_names %>% dplyr::filter(name == .x) %>% dplyr::pull(description)} else {
                                "No summary"}
                                }
                                <br>
                                <a target="_blank" href="?show=gene&query={.x}">Link</a>
                                </center>
                                </p>'))


    # colors used within the network
    geneColor <- color_set_gene_alpha[2]
    proteinColor <- color_set_protein_alpha[2]
    metaboliteColor <- color_set_compound_alpha[2]
    positiveColor <-"rgba(173, 103, 125, 0.8)"
    negativeColor <- "rgba(12, 35, 50, 0.8)"
    connectedColor <- "rgba(255, 255, 255, 0.8)" #formerly purple "rgba(84, 64, 151, 0.8)"
    borderColor <- "rgba(204, 204, 204, 0.8)" #(gray80), formerly white 255, 255, 255
    edgeColor <- "rgba(84, 84, 84, 1)"

    # Physics parameters defining the visNetwork
    iter <- 150 # number of iterations to perform of stablization before display
    gravity <- 0.5
    damping <- 0.11
    timestep <- 0.25

    visNetwork::visIgraph(igraph = metabolic_network) %>%
      visOptions(highlightNearest = list(enabled = T)) %>%
      visGroups(groupname = "Query",
                color = list(background = dplyr::if_else(input$type == "gene", geneColor, metaboliteColor),
                             border =borderColor,
                             highlight = dplyr::if_else(input$type == "gene", geneColor, metaboliteColor),
                             hover = dplyr::if_else(input$type == "gene", geneColor, metaboliteColor)),
                shape = dplyr::if_else(input$type == "gene", 'square', 'dot'),
                borderWidth = 2) %>%
      visGroups(groupname = "Protein",
                color = list(background = proteinColor,
                             border =borderColor,
                             highlight = proteinColor,
                             hover = proteinColor ),
                shape='square',
                borderWidth = 2) %>%
      visGroups(groupname = "Metabolite",
                color = list(background = metaboliteColor,
                             border = borderColor,
                             highlight = metaboliteColor,
                             hover = metaboliteColor),
                shape='dot',
                borderWidth = 2) %>%
      visEdges(color = edgeColor, smooth = F) %>%
      visNodes(scaling = list(min = 10, max =30)) %>%
      visEvents(type = "once",
                startStabilizing = "function() {this.moveTo({scale:0.5})}") %>%
      visPhysics(barnesHut = list(damping = damping, centralGravity = gravity),
                 timestep = timestep,
                 stabilization = FALSE)
  }
  #error handling
  tryCatch(make_bipartite_graph_raw(),
           error = function(e){
             message(e)
             make_empty_graph()
           })
}

