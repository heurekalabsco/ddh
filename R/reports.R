#build a summary table from fun_text.R funs() for report
make_gene_summary <- function(data_values) { #do I need to carry over summary table vars?
  if(data_values$subtype == "gene"){
    summary_table <- tibble(
      identifier = summary_gene(summary_table = gene_summary, data_values, var = "approved_symbol"),
      name = summary_gene(summary_table = gene_summary, data_values, var = "approved_name"),
      summary = summary_gene(summary_table = gene_summary, data_values, var = "entrez_summary"))
  } else if (data_values$subtype == "pathway") {
    summary_table <- tibble(
      identifier = paste0("GO:", summary_pathway(summary_table = pathways, data_values, var = "go")),
      name = summary_pathway(summary_table = pathways, data_values, var = "pathway"),
      summary = summary_pathway(summary_table = pathways, data_values, var = "def"))
  } else { #gene_list
    summary_table <- tibble(
      identifier = c("custom gene list"),
      name = summary_gene_list(summary_table = gene_summary, data_values),
      summary = c("User defined gene input list"))
  }
  return(summary_table)
}

#tests
#make_gene_summary(data_values = list(query="SIRT4", type="gene", subtype = "gene", content=c("SIRT4")))
#make_summary(data_values = list(query="0060148", type="pathway"))
#make_summary(data_values = list(query="SDHA,SDHB", content=c("SDHA", "SDHB"), type="gene_list"))

make_compound_summary <- function(data_values) {
  if(data_values$subtype == "compound"){
    summary_table <- tibble(
      identifier = summary_compound(input = data_values, var = "name"),
      name = summary_compound(input = data_values, var = "moa"),
      summary = summary_compound(input = data_values, var = "description"))
  } else if (data_values$subtype == "moa") {
    summary_table <- tibble(
      identifier = data_values$query,
      name =  c(" (moa)"),
      summary = summary_compound(input = data_values, var = "name"))
  } else { #compound_list
    summary_table <- tibble(
      identifier = c("custom compound list"),
      name = summary_compound_list(input = data_values$content),
      summary = c("User defined gene input list"))
  }
  return(summary_table)
}

make_cell_summary <- function(data_values) {
  if(data_values$type == "cell"){
    summary_table <- tibble(
      identifier = summary_cell(summary_table = expression_names, data_values, var = "cell_line"),
      name = summary_cell(summary_table = expression_names, data_values, var = "cell_line"),
      summary = summary_cell(summary_table = expression_names, data_values, var = "cell_line"))
  } else if (data_values$type == "lineage") {
    summary_table <- tibble(
      identifier = summary_lineage(summary_table = expression_names, data_values, var = "cell_line"),
      name = summary_lineage(summary_table = expression_names, data_values, var = "cell_line"),
      summary = summary_lineage(summary_table = expression_names, data_values, var = "cell_line"))
  } else { #cell_list
    summary_table <- tibble(
      identifier = c("custom cell list"),
      name = summary_cell_list(summary_table = expression_names, data_values),
      summary = c("User defined cell line input list"))
  }
  return(summary_table)
}

#render in temp dir replaces usual render function
render_rmarkdown_in_tempdir <- function(data_values,
                                        rmd_path,
                                        output_file,
                                        envir = parent.frame(),
                                        private_var) {
  # The rmd_path variable must be an absolute path.

  # make sure the base report directory exists
  report_base_dir = here::here("report")
  if (!file.exists(report_base_dir)) {
    dir.create(report_base_dir)
  }
  # determine the filename of the Rmd file we will use for rendering
  rmd_filename <- basename(rmd_path)
  # create a temporary directory and make it our working directory
  temp_dir <- tempfile(pattern="tmpdir", tmpdir=report_base_dir)
  dir.create(temp_dir)
  owd <- setwd(temp_dir)
  on.exit(setwd(owd))
  on.exit(unlink(temp_dir, recursive = TRUE))
  # copy the Rmd file into our temporary(current) directory
  file.copy(rmd_path, rmd_filename, overwrite = TRUE)

  #good file names
  good_file_name <- data_values$query
  if (data_values$subtype == "gene_list" | data_values$subtype == "compound_list" | data_values$type == "cell_list") {
    good_file_name <- paste0("custom_", paste0(data_values$content, collapse="_"))
  }

  #zip
  output_html_filename <- paste0(good_file_name, "_REPORT.html")
  zip_filenames <- c(output_html_filename)

  #CODE TO GET NETWORK INTO REPORT
  if(private_var == TRUE & data_values$type != "cell" & data_values$type != "cell_list") {
    # bring in network from parent environment
    network <- get("network", envir = envir)
    network_filename <- paste0(good_file_name, "_", "graph")

    # Save base network with legend added for the final zip
    visSave(network %>% visLegend(position = "right", width = .25, zoom = F),
            file = paste0(network_filename, "_interactive.html"))
    zip_filenames <- append(zip_filenames, paste0(network_filename, "_interactive.html"))

    # assign filename information to envir so the network can be found within the report_gene.Rmd file
    assign("network_path", paste0(network_filename, "_interactive.html"), envir = envir)
  }

  #PULL IN RENDER()
  rmarkdown::render(rmd_filename,
                    output_file = output_html_filename,
                    envir = envir)

  # get the names of all the items included for rendering
  for (name in names(envir)) {
    env_item = envir[[name]]
    # if the env_item is a plot
    if ("plot_env" %in% names(env_item)) {
      # save the plot to a png
      plot_filename <- paste0(good_file_name, "_", name, ".png")
      # custom heights for lineage plots
      if (name == "lineage") {
        ggsave(plot_filename, width = 12, height = 10.5, env_item, dpi = 300, type = "cairo")
      }
      if (name == "sublineage") {
        ggsave(plot_filename, env_item, width = 12, height = 22, dpi = 300, type = "cairo")
      }
      # dynamic height depending on # of genes for cellbins plot
      if (name == "cellbins") {
        ggsave(plot_filename, env_item, width = 12, dpi = 300, type = "cairo",
               height = length(unique(data_values$content)) + 3, limitsize = FALSE)
      }
      # smaller plots for anatogram and network
      if (name == "cellanatogram") {
        ggsave(plot_filename, env_item, width = 8, height = 7, dpi = 300, type = "cairo")
      }
      if (name == "cellanatogram_facet") {
        n <- length(unique(data_values$content))
        h <- (n + 2) %/% 3 * 4
        ## if less than 3 genes the plot should still span the full page
        ## so we need to adjust the height accordingly
        if(n < 3) h <- n * 2
        ggsave(plot_filename, env_item, dpi = 300, type = "cairo",
               width = 12, height = h, limitsize = FALSE)
      }
      if (name == "graph") {
        ggsave(plot_filename, env_item, width = 12, height = 10.5, dpi = 300, type = "cairo")
      }
      ## scale protein size plot based on number of genes
      if (name == "proteinsize") {
        ggsave(plot_filename, env_item, width = 10,
               height = length(unique(data_values$content)) * 1.36 + .9,
               dpi = 300, type = "cairo")
      }
      ## scale expression plots based on number of genes
      if (name %in% c("expression_gene_plot", "expression_protein_plot")) {
        ggsave(plot_filename, env_item, width = 10,
               height = length(unique(data_values$content)) * .7 + .8,
               dpi = 300, type = "cairo")
      }
      # landscape aspect ratio for all other plots
      if (!name %in% c("lineage", "sublineage", "cellbins", "cellanatogram",
                       "cellanatogram_facet", "graph", "proteinsize",
                       "expression_gene_plot", "expression_protein_plot")) {
        ggsave(plot_filename, env_item, width = 12, height = 7.5, dpi = 300, type = "cairo")
      }
      #include prerendered images
      structure_filepath <- make_structure(input = data_values)
      if(!is.null(structure_filepath)){ #error checking
        structure_filename <- basename(structure_filepath)

        file.copy(from = structure_filepath,
                  to = structure_filename, overwrite = TRUE)

        zip_filenames <- append(zip_filenames, structure_filename)
      }

      # include the plot png in the zip download
      zip_filenames <- append(zip_filenames, plot_filename)
    }
  }
  zip(zipfile = output_file, files = zip_filenames)
}

#specific instructions to render reports based on query type and report template
render_gene_report <- function(data_values, output_file, privateMode) {
  #no need to declare type here b/c that's done in render_report_to_file() below
  num <- dplyr::n_distinct(achilles_long$X1)
  summary <- make_gene_summary(data_values)
  ideogram <- make_ideogram(input = data_values)
  proteinsize_plot <- make_proteinsize(protein_data = proteins, input = data_values, card = FALSE)
  structure_path <- make_structure(input = data_values)
  pubmed_plot <- make_pubmed(input = data_values)
  pubmed_table <- make_pubmed_table(input = data_values)
  cellanatogram <- make_cellanatogram(cellanatogram_data = subcell, input = data_values)
  cellanatogram_facet <- make_cellanatogramfacet(cellanatogram_data = subcell, input = data_values)
  cellanatogram_table <- make_cellanatogram_table(cellanatogram_data = subcell, gene_symbol = data_values$content)
  expression_gene_plot <- make_cellexpression(input = data_values, var = "gene")
  expression_gene_table <- make_expression_table(input = data_values, var = "gene")
  expression_protein_plot <- make_cellexpression(input = data_values, var = "protein")
  expression_protein_table <- make_expression_table(input = data_values, var = "protein")
  expression_cellgeneprotein <- make_cellgeneprotein(input = data_values)
  colorful_male <- make_male_anatogram(input = data_values)
  colorful_female <- make_female_anatogram(input = data_values)
  tissue_plot <- make_tissue(tissue_data = tissue, input = data_values)
  humananatogram_table <- make_humananatogram_table(humananatogram_data = tissue, input = data_values)
  celldeps <- make_celldeps(input = data_values)
  cellbars <- make_cellbar(input = data_values)
  cellbins <- make_cellbins(input = data_values)
  genecorrelation <- make_correlation(input = data_values)
  lineage <- make_lineage(input = data_values)
  sublineage <- make_sublineage(input = data_values)
  dep_table <- make_dep_table(input = data_values)
  dep_top <- make_top_table(input = data_values)
  flat_top_complete <- make_enrichment_top(input = data_values)
  dep_bottom <- make_bottom_table(input = data_values)
  flat_bottom_complete <- make_enrichment_bottom(input = data_values)
  expdep <- make_expdep(input = data_values)
  drug_cor_table <- make_gene_drugs_cor_table(gene_symbol = data_values$content)
  if (privateMode == TRUE) {
    #make some private data
    network <- make_graph(input = data_values, threshold = 10, deg = 2, corrType = "Positive", displayHeight = '80vh', displayWidth = '100%')
    #signature
    radial_aa_plot <- make_radial(input = data_values, card = FALSE)
    radial_aa_plot_clust <- make_radial(input = data_values, cluster = TRUE, card = FALSE)
    clustering_table <- make_clustering_table(gene_symbol = data_values$content)
    clustering_enrichment_table_bp <- make_clustering_enrichment_table(gene_symbol = data_values$content, ontology = "BP")
    clustering_enrichment_table_mf <- make_clustering_enrichment_table(gene_symbol = data_values$content, ontology = "MF")
    clustering_enrichment_table_cc <- make_clustering_enrichment_table(gene_symbol = data_values$content, ontology = "CC")
    cluster_enrich_plot_bp <- make_cluster_enrich(input = data_values, ontology = "BP")
    cluster_enrich_plot_mf <- make_cluster_enrich(input = data_values, ontology = "MF")
    cluster_enrich_plot_cc <- make_cluster_enrich(input = data_values, ontology = "CC")
    #drug_table <- tryCatch(make_gene_drugs_table(gene_symbol = data_values$content), error = function(x){as.data.frame("No data available")})
    drug_table <- make_gene_drugs_table(gene_symbol = data_values$content)
    metabolite_table <- make_metabolite_table(input = data_values)
    bipartite_graph <- make_bipartite_graph(input = data_values)
  }
  #this calls the rmd
  render_rmarkdown_in_tempdir(data_values,
                              here::here("code", "report_gene.Rmd"),
                              output_file,
                              private_var = privateMode)
}

# render_gene_report(data_values = list(query="HDDC3", type="gene", subtype="gene", content=c("HDDC3")), output_file=here::here("/report/HDDC3.zip"), privateMode = TRUE)
# render_gene_report(data_values = list(query="ROCK1", type="gene", subtype="gene", content=c("ROCK1")), output_file="ROCK1_manual.zip",  privateMode = TRUE)
# #render_gene_report(input = "0060148", type = "pathway", output_file = "0060148.zip")
# #render_gene_report(input = c("GSS", "SST"), type = "gene_list")

render_compound_report <- function(data_values, output_file, privateMode) {
  #no need to declare type here b/c that's done in render_report_to_file() below
  num <- dplyr::n_distinct(prism_names$name)
  summary <- make_compound_summary(data_values)
  pubmed_plot <- make_pubmed(input = data_values)
  pubmed_table <- make_pubmed_table(input = data_values)
  #strucutre <- make_molecule_structure(input = data_values)
  if (privateMode == TRUE) {
    #some private data
    celldeps <- make_celldeps(input = data_values)
    cellbins <- make_cellbins(input = data_values)
    lineage <- make_lineage(input = data_values)
    sublineage <- make_sublineage(input = data_values)
    dep_table <- make_dep_table(input = data_values)
    dep_top <- make_compound_table(input = data_values, top = TRUE)
    dep_bottom <-make_compound_table(input = data_values, top = FALSE)
    network <- make_graph(input = data_values, threshold = 10, deg = 2, corrType = "Positive", displayHeight = '80vh', displayWidth = '100%')
  }
  #this calls the rmd
  render_rmarkdown_in_tempdir(data_values,
                              here::here("code", "report_compound.Rmd"),
                              output_file,
                              private_var = privateMode)
}
# render_compound_report(data_values = list(query="aspirin", type="compound", compound=c("aspirin")), output_file=here::here("/report/dataTest.zip"))

render_cell_report <- function(data_values, output_file, privateMode) {
  #no need to declare type here b/c that's done in render_report_to_file() below
  num <- dplyr::n_distinct(achilles_long$X1)
  summary <- make_cell_summary(data_values)
  if (privateMode == TRUE) {
    celldeps <- make_celldeps(input = data_values)
    cellbins <- make_cellbins(input = data_values)
    expression_gene_plot <- make_cellexpression(input = data_values, var = "Gene")
    expression_gene_table <- make_expression_table(input = data_values, var = "Gene")
    expression_protein_plot <- make_cellexpression(input = data_values, var = "Protein")
    expression_protein_table <- make_expression_table(input = data_values, var = "Protein")
    expression_cellgeneprotein <- make_cellgeneprotein(input = data_values)
    dep_table <- make_dep_table(input = data_values)
  }
  #this calls the rmd
  render_rmarkdown_in_tempdir(data_values,
                              here::here("code", "report_cell_line.Rmd"),
                              output_file,
                              private_var = privateMode)
}

#logic to matching query type to rendered content and set privateMode
render_report_to_file <- function(data_values,
                                  file,
                                  privateMode) {
  if (data_values$type == "gene") {
    if (privateMode == TRUE) {
      render_gene_report(data_values, output_file = file, privateMode = TRUE)
    } else {
      render_gene_report(data_values, output_file = file, privateMode = FALSE)
    }
  } else if (data_values$type == "compound") {
    if (privateMode == TRUE) {
      render_compound_report(data_values, output_file = file, privateMode = TRUE)
    } else {
      render_compound_report(data_values, output_file = file, privateMode = FALSE)
    }
  } else if (data_values$type == "cell") {
    if (privateMode == TRUE) {
      render_cell_report(data_values, output_file = file, privateMode = TRUE)
    } else {
      render_cell_report(data_values, output_file = file, privateMode = FALSE)
    }
  } else {
    stop("no report for you!")
  }
}

#render_report_to_file(data_values = list(query="ROCK1", type="gene", subtype="gene", content=c("ROCK1")), file=here::here("/report/ROCK1_manual_full.zip"),  privateMode = TRUE)
