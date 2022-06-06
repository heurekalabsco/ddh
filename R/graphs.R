## EMPTY GRAPH -----
#' Empty Graph Graph
#'
#' \code{make_empty_graph} returns an image of ...
#'
#' This is a graph function that takes a gene name and returns a empty graph graph
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a empty graph graph. If an error is thrown, then will return an empty graph.
#'
#' @examples
#' make_empty_graph()
#' \dontrun{
#' make_empty_graph()
#' }
make_empty_graph <- function(type = "gene") {
  if(type == "gene") {
    queryColor <- color_set_gene_alpha[2]
  } else if(type == "compound") {
    queryColor <- color_set_compound_alpha[2]
  } else {
    stop("declare your type")
  }
  #set params, copied from below
  borderColor <- "rgba(204, 204, 204, 0.8)" #(gray80), formerly white 255, 255, 255
  displayHeight = '90vh'
  displayWidth = '100%'
  #make empty nodes table
  nodes_empty = tibble(id = 0,
                       name = "Empty",
                       group = "Query",
                       title = "Not enough data to build a network")
  visNetwork(nodes = nodes_empty,
             width = displayWidth,
             height = displayHeight) %>%
    visGroups(groupname = "Query",
              color = list(background = queryColor, border = borderColor, highlight = queryColor, hover = queryColor),
              shape='dot',
              borderWidth = 2)
}

## SETUP GRAPH ----------------------------------------------------------------------
setup_graph <- function(toptable_data = master_top_table,
                        bottomtable_data = master_bottom_table,
                        compound_data = prism_cor_nest,
                        input_list = list(), #changed name here to prevent var naming overlap for nested funs()
                        setup_threshold,
                        setup_corrType) {

  if(input_list$type == "gene") {
    #generate data
    #either find top/bottom correlated genes if given single gene, or take list to fill gene_list
    if(length(input_list$content) == 1){
      #find top and bottom correlations for fav_gene
      query <- input_list$content
      if(setup_corrType == "Positive" | setup_corrType == "Positive and Negative") {
        top <-
          toptable_data %>%
          dplyr::filter(fav_gene %in% input_list$content) %>%
          tidyr::unnest(data) %>%
          dplyr::arrange(dplyr::desc(r2)) %>%
          dplyr::slice(1:setup_threshold) %>%
          dplyr::pull("gene")
      }
      if(setup_corrType == "Negative" | setup_corrType == "Positive and Negative") {
        bottom <-
          bottomtable_data %>%
          dplyr::filter(fav_gene %in% input_list$content) %>%
          tidyr::unnest(data) %>%
          dplyr::arrange(r2) %>%
          dplyr::slice(1:setup_threshold)%>%
          dplyr::pull("gene")
      }
    } else {
      query <- input_list$content
      top <- input_list$content #set to query here, to pull top correlated genes, reset below
      bottom <- input_list$content  #set to query here, to pull bottom correlated genes, reset below
    }
  } else if(input_list$type == "compound") {
    #do this
    if(length(input_list$content) == 1){
      #find top and bottom correlations for fav_drug
      query <- input$content
      if(setup_corrType == "Positive" | setup_corrType == "Positive and Negative") {
        top <-
          compound_data %>%
          dplyr::filter(fav_drug %in% input_list$content) %>%
          tidyr::unnest("data") %>%
          dplyr::ungroup() %>%
          dplyr::filter(r2 > prism_cor_upper) %>% #filter(., r2 < prism_cor_lower)
          dplyr::arrange(desc(r2)) %>%
          dplyr::slice(1:setup_threshold) %>%
          dplyr::pull("name")
      }
      if(setup_corrType == "Negative" | setup_corrType == "Positive and Negative") {
        bottom <-
          compound_data %>%
          dplyr::filter(fav_drug %in% input_list$content) %>%
          tidyr::unnest("data") %>%
          dplyr::ungroup() %>%
          dplyr::filter(r2 < prism_cor_lower) %>%
          dplyr::arrange(r2) %>%
          dplyr::slice(1:setup_threshold) %>%
          dplyr::pull("name")
      }
    } else {
      query <- input_list$content
      top <- input_list$content
      bottom <- input_list$content
    }
  } else {
    stop("declare your type")
  }

  #table maker function
  make_graph_table <- function(fun_input_list = list(),
                               content,
                               setup_threshold,
                               top = TRUE) {
    if(fun_input_list$type == "gene") {
      filter_var <- rlang::sym("fav_gene")
      rename_var <- rlang::sym("gene")
      if(top == TRUE) {
        message_var <- "top"
        table_var <- toptable_data
        origin_var <- "pos"
      } else {
        message_var <- "bottom"
        table_var <- bottomtable_data
        origin_var <- "neg"
      }
    } else if(fun_input_list$type == "compound") {
      filter_var <- rlang::sym("fav_drug")
      rename_var <- rlang::sym("name")
      if(top == TRUE) {
        message_var <- "top"
        table_var <- compound_data
        origin_var <- "pos"
      } else {
        message_var <- "bottom"
        table_var <- compound_data
        origin_var <- "neg"
      }
    } else {
      stop("declare your type")
    }

    message(glue::glue('Getting {message_var} correlations from {content}'))

    if(content %in% table_var[[1]]){ #check to see if gene query is in table (either fav_gene or fav_drug)
      related_table <-
        table_var %>%
        dplyr::filter(!!filter_var == content) %>%
        tidyr::unnest(data) %>%
        dplyr::ungroup(.) %>%
        {if (top == TRUE) dplyr::arrange(., desc(r2)) else dplyr::arrange(., r2)} %>%
        dplyr::slice(1:setup_threshold) %>%
        dplyr::mutate(x = content,
                      origin = origin_var) %>%
        dplyr::rename(y = !!rename_var) %>%
        dplyr::select(x, y, r2, origin)
      return(related_table)}
  }

  #make empty tibble
  dep_network <- tibble()

  # make the correct graph including only correlations of the designated type
  if(setup_corrType == "Positive"){
    network_list <- unique(top)
    bottom <- NULL #reset unused var to NULL so factors and labels work
    #this takes the genes from the top, and pulls them to feed them into a for loop
    for (i in network_list){
      dep_top_related <- make_graph_table(fun_input_list = input_list,
                                          content = i,
                                          setup_threshold,
                                          top = TRUE)
      #each temp object is bound together
      dep_network <-
        dep_network %>%
        bind_rows(dep_top_related)
    }
  } else if(setup_corrType == "Positive and Negative"){
    network_list <- unique(c(top, bottom))
    for (i in network_list){
      dep_top_related <- make_graph_table(fun_input_list = input_list,
                                          content = i,
                                          setup_threshold,
                                          top = TRUE)
      dep_bottom_related <- make_graph_table(fun_input_list = input_list,
                                             content = i,
                                             setup_threshold,
                                             top = FALSE)

      #each temp object is bound together, and then bound to the final df for graphing
      dep_related <-
        dep_top_related %>%
        bind_rows(dep_bottom_related)

      dep_network <-
        dep_network %>%
        bind_rows(dep_related)
    }
  } else if(setup_corrType == "Negative"){
    network_list <- unique(bottom)
    top <- NULL #reset unused var to NULL so factors and labels work
    for (i in network_list){
      dep_bottom_related <- make_graph_table(fun_input_list = input,
                                             content = i,
                                             setup_threshold,
                                             top = FALSE)

      #each temp object is bound together
      dep_network <-
        dep_network %>%
        bind_rows(dep_bottom_related)
    }
  } else {
    stop("delcare your corrType")
  }
  if(length(input_list$content) > 1) { #this reassigns top and bottom ids from multi-gene queries
    query <- input_list$content
    top <- dep_network %>% filter(origin == "pos") %>% pull(y)
    bottom <- dep_network %>% filter(origin == "neg") %>% pull(y)
  }
  #fun now returns a list, to preserve some info that gets passed on to rest make_graph()
  return(list(type = setup_corrType,
              query_id = query,
              top_id = top,
              bottom_id = bottom,
              network_l = network_list,
              df = dep_network))
}

#tests
#tmp <- setup_graph(input_list = list(type = "gene", content = "ROCK1"), setup_corrType = "Positive", setup_threshold = 10)
#tmp1 <- setup_graph(input_list = list(type = "gene", content = "ROCK1"), setup_corrType = "Positive and Negative", setup_threshold = 10)
#tmp2 <- setup_graph(input_list = list(type = "gene", content = c("ROCK1", "ROCK2")), setup_corrType = "Positive", setup_threshold = 10)
#tmp4 <- setup_graph(input_list = list(type = "compound", content = "ADP"), setup_corrType = "Positive and Negative", setup_threshold = 10)

#' Graph Graph
#'
#' \code{make_graph} returns an image of ...
#'
#' This is a graph function that takes a gene name and returns a graph graph
#'
#' @param input Expecting a list containing type and content variable.
#'


## MAKE GRAPH ----------------------------------------------------------------------
#' Create network graph visualization using visNetwork
#'
#' This function takes in dependency correlations and a gene query list to then output a dependency network graph
#' visualization containing the top/bottom threshold for each of the top/bottom threshold of the gene query list
#' using visNetwork.
#'
#' @param toptable_data A tibble of genes and their associated top correlated genes
#' @param bottomtable_data A tibble of genes and their associated bottom correlated genes
#' @param input A list containing character vector of gene_symbols used to create network graph
#' @param threshold A numerical representing the number of genes to pull from top and bottom tables
#' @param deg A numerical representing the minimum number of connections for a gene to be connected to the network
#' @param corrType A string that describes what type of correlations to include, options are: "Positive and Negative", "Positive", or "Negative"
#' @param displayHeight Default to "90vh". The height of the network in pixels("500px"), as a percentage("100%"), or as a percentage of the viewport("70vh", where 70 represents 70% of the viewport)
#' @param displayWidth Default to "100%". The width of the network in pixels("500px"), as a percentage("100%"), or as a percentage of the viewport("70vh", where 70 represents 70% of the viewport)
#' @param tooltipLink Boolean to denote whether or not to include a link in the tooltip for a gene. Default to false.
#'
#' @return Outputs a complete network graph. If an error is thrown, then will return an empty graph.
#' @export
#' @examples
#' make_graph(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_graph(input = list(type = "gene", content = "ROCK1"))
#' make_graph(input = list(type = "gene", content = "ROCK1"), card = TRUE)
#' make_graph(input = list(type = "gene", content = "ROCK1"), corrType = "Positive")
#' make_graph(input = list(type = "pathway", query = "1902965", content = c("RDX", "ROCK2", "DTX3L", "MSN", "SORL1", "EZR")), corrType = "Positive")
#' make_graph(input = list(type = "gene", content = c("ROCK1", "ROCK2")))
#' make_graph(input = list(type = "gene", content = "DTX3L"), corrType = "Negative") # disconnected query gene
#' make_graph(input = list(type = "compound", content = "aspirin"), corrType = "Negative")

#' \dontrun{
#' make_graph(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_graph <- function(toptable_data = master_top_table,
                       bottomtable_data = master_bottom_table,
                       compound_data = prism_cor_nest,
                       input = list(),
                       threshold = 10,
                       deg = 2,
                       corrType = "Positive",
                       displayHeight = '90vh',
                       displayWidth = '100%',
                       tooltipLink = FALSE,
                       card = FALSE) {
  make_graph_raw <- function() {
    if(input$type == "gene") {
      queryColor <- color_set_gene_alpha[2]
    } else if(input$type == "compound") {
      queryColor <- color_set_compound_alpha[2]
    } else {
      stop("declare your type")
    }

    #get dep_network object
    dep_network_list <- setup_graph(input_list = input,
                                    setup_corrType = corrType,
                                    setup_threshold = threshold)
    #dep_network_list <<- dep_network_list #for testing, to see what I'm getting back out

    #add check for no genes list
    #if(nrow(dep_network_list$df) == 0) {return("Graph cannot be built")}
    #comment out so throws error with tryCatch()

    if(length(input$content) == 1){
      if(corrType == "Positive") {
        group_var <- c("Query", "Positive", "Connected")
      } else if(corrType == "Positive and Negative") {
        group_var <- c("Query", "Positive", "Negative", "Connected")
      } else if(corrType == "Negative") {
        group_var <- c("Query", "Negative", "Connected")
      } else {
        stop("declare your corrType")}
    } else {
      if(corrType == "Positive") {
        group_var <- c("Query", "Positive")
      } else if(corrType == "Positive and Negative") {
        group_var <- c("Query", "Positive", "Negative")
      } else if(corrType == "Negative") {
        group_var <- c("Query", "Negative")
      } else {
        stop("declare your corrType")}
    }

    #make graph
    graph_network <-
      tidygraph::as_tbl_graph(dep_network_list$df)

    nodes <-
      as_tibble(graph_network) %>%
      rowid_to_column("id") %>%
      mutate(degree = igraph::degree(graph_network),
             group = case_when(name %in% dep_network_list$query_id == TRUE ~ "Query", #could use input$content or input$content
                               name %in% dep_network_list$top_id == TRUE ~ "Positive",
                               name %in% dep_network_list$bottom_id == TRUE ~ "Negative",
                               TRUE ~ "Connected"),
             group = as_factor(group),
             group = fct_relevel(group, group_var))  %>%
      arrange(group)

    links <- graph_network %>%
      activate(edges) %>% # %E>%
      as_tibble()

    # determine the nodes that have at least the minimum degree
    nodes_filtered <-
      nodes %>%
      filter(degree >= deg) %>%  #input$degree
      as.data.frame

    # filter the edge list to contain only links to or from the nodes that have the minimum or more degree
    links_filtered <-
      links %>%
      filter(to %in% nodes_filtered$id & from %in% nodes_filtered$id) %>%
      as.data.frame

    # re-adjust the from and to values to reflect the new positions of nodes in the filtered nodes list
    links_filtered$from <- match(links_filtered$from, nodes_filtered$id) - 1
    links_filtered$to <- match(links_filtered$to, nodes_filtered$id) - 1

    #check to see if setting degree removed all links; if so, then throws error, so this fills a dummy links_filtered df to plot only nodes
    if(nrow(links_filtered) == 0) {
      if(corrType == "Negative"){
        links_filtered <- tibble("from" = -1, "to" = -1, "r2" = 1, "origin" = "neg")
      }    else{
        links_filtered <- tibble("from" = -1, "to" = -1, "r2" = 1, "origin" = "pos")}
    }

    # shift id values properly to work with visNetwork
    nodes_filtered <-
      nodes_filtered %>%
      dplyr::mutate(id=0:(dim(nodes_filtered)[1]-1))

    # recreate the network using filtered edges to use degree function
    graph_network_filtered <-
      igraph::graph_from_data_frame(links_filtered, directed = F) %>%
      igraph::simplify()
    degVec <- igraph::degree(graph_network_filtered)

    # make node size a function of their degree in the current network
    if(names(degVec)[1] == "-1"){ # if there are no links remaining then assign degree of each node to 0
      nodes_filtered <- nodes_filtered %>%
        mutate(value = 0)
    } else if(length(input$content > 1) & length(degVec) != dim(nodes_filtered)[1]){ # handle cases where some query genes are disconnected
      genesWithConnections <- degVec %>%
        names() %>%
        as.numeric()
      nodes_filtered <- nodes_filtered %>%
        add_column(value = 0)
      for(gene in 1:dim(nodes_filtered)[1]){
        if(nodes_filtered[gene,"id"] %in% genesWithConnections){
          nodes_filtered[gene,"value"] <- degVec[nodes_filtered[gene,"id"] %>% toString()]
        }
      }
    }else{
      nodes_filtered <- nodes_filtered[degVec %>% names() %>% as.numeric() +1,] %>%
        mutate(value = degVec)
    }

    # make sure query gene is at the start so legend shows proper order
    nodes_filtered <- nodes_filtered %>%
      arrange(id)

    #check to see if query gene is missing; if so, then adds a dummy so it shows up on graph, but disconnected
    disconnected <- F
    #could simplify this if/else with common input name, b/c only here to get input$content vs. input$content
    if(input$type == "gene") {
      if(sum(str_detect(nodes_filtered$group, "Query")) == 0){
        dummy <- tibble("id" = max(nodes_filtered$id) + 1, "name" = input$content, "degree" = 1, "group" = "Query", "value" = 0)
        nodes_filtered <- bind_rows(dummy, nodes_filtered)
        disconnected <- T
      }
    } else if(input$type == "compound") {
      if(sum(str_detect(nodes_filtered$group, "Query")) == 0){
        dummy <- tibble("id" = max(nodes_filtered$id) + 1, "name" = input$content, "degree" = 1, "group" = "Query", "value" = 0)
        nodes_filtered <- bind_rows(dummy, nodes_filtered)
        disconnected <- T
      }
    } else {
      stop("declare your type")
    }

    # get the approved_name of each gene from the gene_summary table - will be added to nodes tibble for tooltip
    nameTable <- tibble(name=character())

    if(input$type == "gene") {
      for(gene in nodes_filtered$name){
        newVal <- gene_summary %>%
          dplyr::filter(approved_symbol==gene) %>%
          dplyr::pull(approved_name)
        if(length(newVal)==0){
          nameTable <- add_row(nameTable, name = "No gene summary available")# handles cases where the gene is not in the gene summary table
        } else{
          nameTable <- add_row(nameTable, name=newVal)
        }
      }
    } else if(input$type == "compound") {
      for(drug in nodes_filtered$name){
        newVal <- prism_meta %>%
          dplyr::filter(name==drug) %>%
          dplyr::pull(moa)
        if(length(newVal)==0){
          nameTable <- add_row(nameTable, name = "No drug MOA available")# handles cases where the gene is not in the gene summary table
        } else{
          nameTable <- add_row(nameTable, name=newVal)
        }
      }
    } else {
      stop("declare your type")
    }

    # add title information (tooltip that appears on hover)
    if(!tooltipLink){ # Do not form a url when just making the standalone graph while testing or making reports
      nodes_filtered <- nodes_filtered %>%
        dplyr::mutate(title=paste0("<center><p>", nodes_filtered$name,"<br>",nameTable$name, '</p>'),
                      label = nodes_filtered$name )
    }else{
      if(input$type == "gene") {
        nodes_filtered <- nodes_filtered %>%
          dplyr::mutate(title=paste0("<center><p>", nodes_filtered$name,"<br>",nameTable$name ,'<br><a target="_blank" href="?show=gene&query=',nodes_filtered$name,'">Gene Link</a></p>'),
                        label = nodes_filtered$name )
      } else if(input$type == "compound") {
        nodes_filtered <- nodes_filtered %>%
          dplyr::mutate(title=paste0("<center><p>", nodes_filtered$name,"<br>",nameTable$name ,'<br><a target="_blank" href="?show=compound&query=',nodes_filtered$name,'">Drug Link</a></p>'),
                        label = nodes_filtered$name )
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
    if(corrType == "Positive and Negative"){
      visNetwork(nodes = nodes_filtered, edges = links_filtered, width = displayWidth, height = displayHeight) %>%
        visOptions(highlightNearest = list(enabled = T)) %>%
        visGroups(groupname = "Query", color = list(background = queryColor, border =borderColor, highlight = queryColor, hover = queryColor ), shape='dot', borderWidth = 2) %>%
        visGroups(groupname = "Positive", color = list(background = positiveColor, border = borderColor, highlight = positiveColor, hover = positiveColor), shape='dot', borderWidth = 2) %>%
        visGroups(groupname = "Negative", color = list(background = negativeColor, border = borderColor, highlight = negativeColor, hover = negativeColor), shape='dot', borderWidth = 2) %>%
        visGroups(groupname = "Connected", color = list(background = connectedColor, border = borderColor, highlight = connectedColor, hover = connectedColor), shape='dot', borderWidth = 2) %>%
        visEdges(color = edgeColor, smooth = F) %>%
        visNodes(scaling = list(min = 10, max =20)) %>%
        visPhysics(barnesHut = list(damping = damping, centralGravity = gravity), timestep = timestep, stabilization = list(iterations = iter)) %>%
        visEvents(stabilizationIterationsDone = stabilizationZoomFn)
      # visConfigure(enabled=TRUE) # use to test out new features to add from visNetwork
    }  else if(corrType == "Positive"){
      visNetwork(nodes = nodes_filtered, edges = links_filtered, width = displayWidth, height = displayHeight) %>%
        visOptions(highlightNearest = list(enabled = T)) %>%
        visGroups(groupname = "Query", color = list(background = queryColor, border =borderColor, highlight = queryColor, hover = queryColor ), shape='dot', borderWidth = 2) %>%
        visGroups(groupname = "Positive", color = list(background = positiveColor, border = borderColor, highlight = positiveColor, hover = positiveColor), shape='dot', borderWidth = 2) %>%
        visGroups(groupname = "Connected", color = list(background = connectedColor, border = borderColor, highlight = connectedColor, hover = connectedColor), shape='dot', borderWidth = 2) %>%
        visEdges(color = edgeColor, smooth = F) %>%
        visNodes(scaling = list(min = 10, max =20)) %>%
        visPhysics(barnesHut = list(damping = damping, centralGravity = gravity), timestep = timestep, stabilization = list(iterations = iter)) %>%
        visEvents(stabilizationIterationsDone = stabilizationZoomFn)
    }  else if(corrType == "Negative"){
      visNetwork(nodes = nodes_filtered, edges = links_filtered, width = displayWidth, height = displayHeight) %>%
        visOptions(highlightNearest = list(enabled = T)) %>%
        visGroups(groupname = "Query", color = list(background = queryColor, border =borderColor, highlight = queryColor, hover = queryColor ), shape='dot', borderWidth = 2) %>%
        visGroups(groupname = "Negative", color = list(background = negativeColor, border = borderColor, highlight = negativeColor, hover = negativeColor), shape='dot', borderWidth = 2) %>%
        visGroups(groupname = "Connected", color = list(background = connectedColor, border = borderColor, highlight = connectedColor, hover = connectedColor), shape='dot', borderWidth = 2) %>%
        visEdges(color = edgeColor, smooth = F) %>%
        visNodes(scaling = list(min = 10, max = 20)) %>%
        visPhysics(barnesHut = list(damping = damping, centralGravity = gravity), timestep = timestep, stabilization = list(iterations = iter)) %>%
        visEvents(stabilizationIterationsDone = stabilizationZoomFn)
    }
  }
  #error handling
  tryCatch(make_graph_raw(),
           error = function(x){
             make_empty_graph()
           })
}

#figure legend
graph_title <- "Network Graph."
graph_legend <- "Each point represents a single gene taken from the top associated genes with the query gene. Genes with only one connection were removed."
graph_legend_list <- "Each point represents one of the queried genes, and then the top and bottom associated genes with it. Genes with only one connection were removed."

#' Bipartite Graph Graph
#'
#' \code{make_bipartite_graph} returns an image of ...
#'
#' This is a graph function that takes a gene name and returns a bipartite graph graph
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a bipartite graph graph. If an error is thrown, then will return an empty graph.
#'
#' @export
#' @examples
#' make_bipartite_graph(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' \dontrun{
#' make_bipartite_graph(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_bipartite_graph <- function(toptable_data = master_top_table,
                                 bottomtable_data = master_bottom_table,
                                 protein_data = hmdb_proteins,
                                 metabolite_data = hmdb_metabolites,
                                 input = list(),
                                 censor = character(), #removes most common metabolites
                                 collapsed = TRUE,
                                 threshold = 10,
                                 corrType = "Positive") {
  make_bipartite_graph_raw <- function() {
    if(input$type == "gene") {
      #get dep_network object
      dep_network_list <- setup_graph(input_list = input,
                                      setup_corrType = corrType,
                                      setup_threshold = threshold)
      #get gene_names
      if(length(input$content) == 1){
        if(corrType == "Positive") {
          genes <- c(dep_network_list$query_id, dep_network_list$top_id)
        } else if(corrType == "Positive and Negative") {
          genes <- c(dep_network_list$query_id, dep_network_list$top_id, dep_network_list$bottom_id)
        } else if(corrType == "Negative") {
          genes <- c(dep_network_list$query_id, dep_network_list$bottom_id)
        } else {
          stop("declare your corrType")}
      } else {
        genes <- c(dep_network_list$query_id)
      }

      #collapsed var
      if(collapsed == TRUE) {collapsed_var <- rlang::sym("data_collapsed")
      } else {collapsed_var <- rlang::sym("data_original")}

      #first two cols to generate graph
      hmdb_network <-
        protein_data %>%
        dplyr::filter(fav_gene %in% genes)

      #check for nrow
      if(nrow(hmdb_network) == 0) {
        return("No associated metabolites")
      } else {
        hmdb_network <-
          hmdb_network %>%
          tidyr::unnest(cols = c(!!collapsed_var)) %>%
          dplyr::ungroup() %>%
          dplyr::select(gene_name, metabolite_name)
      }

      #get genes that have >1 connection for filtering
      connected_genes <-
        hmdb_network %>%
        dplyr::count(gene_name, sort = TRUE) %>%
        dplyr::filter(n >1) %>%
        dplyr::pull(gene_name)

      hmdb_filtered <-
        hmdb_network %>%
        dplyr::filter(gene_name %in% connected_genes)

      #hackish way to reorder legend, which appears based on first instance of  gene in group assignment
      hmdb_query <- #get query genes
        hmdb_filtered %>%
        dplyr::filter(gene_name == dep_network_list$query_id)

      hmdb_notquery <- #get all others
        hmdb_filtered %>%
        dplyr::filter(gene_name != dep_network_list$query_id)

      hmdb_filtered <- #bind rows with query first
        hmdb_query %>%
        dplyr::bind_rows(hmdb_notquery)

      #remove common metabolites; receives character vec from shiny with pre-defined choices
      if(length(censor) > 0) {
        hmdb_filtered <-
          hmdb_filtered %>%
          dplyr::filter(!metabolite_name %in% censor)
      }
    } else if (input$type == "compound") {
      #first two cols to generate graph
      hmdb_network <-
        metabolite_data %>%
        dplyr::filter(fav_metabolite %in% input$content)

      #check for nrow
      if(nrow(hmdb_network) == 0) {
        return("No associated proteins")
      } else {
        hmdb_filtered <-
          hmdb_network %>%
          tidyr::unnest(cols = c(data)) %>%
          dplyr::ungroup() %>%
          dplyr::select(gene_name, metabolite_name) %>%
          dplyr::arrange(gene_name)
      }
    } else {
      stop("delcare your type")
    }
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

    # add title information (tooltip that appears on hover)
    igraph::V(metabolic_network)$title <- igraph::V(metabolic_network)$name %>%
      purrr::map_chr(~ glue::glue('<center>
                                <p>
                                {.x}
                                <br>
                                {if(length(gene_summary %>% filter(approved_symbol == .x) %>% pull(approved_name)) != 0) {
                                gene_summary %>% filter(approved_symbol == .x) %>% pull(approved_name)} else {
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
           error = function(x){"Graph cannot be built"})
}
# Test Cases
# make_bipartite_graph(input = list(type = "gene", content = "ROCK1"),
#                      collapsed = TRUE,
#                      threshold = 10,
#                      corrType = "Positive")
#
# make_bipartite_graph(input = list(type = "gene", content = c("ROCK1", "ROCK2")),
#                      collapsed = TRUE,
#                      threshold = 10,
#                      corrType = "Positive",
#                      censor = c("ADP", "Adenosine triphosphate"))
#
# make_bipartite_graph(input = list(type = "compound", content = "adp"))
