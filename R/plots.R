## BARCODE PLOT --------------------------------------------------------------------
#' Barcode Plot
#'
#' The barcode image shows the location of the query gene(s) on human chromosomes in the style of a barcode.
#'
#' @param input Expecting a list containing a content variable.
#' @param card A boolean that sets which image should be pulled for display
#' @return If no error, then returns a barcode image url.
#'
#' @export
#' @examples
#' make_barcode(input = list(type = "gene", query = "ROCK1", content = "ROCK1"))
#' make_barcode(input = list(content = "ROCK1"))
#' \dontrun{
#' make_barcode(input = list(type = "gene", content = "ROCK1"), card = TRUE)
#' }
make_barcode <- function(input = list(),
                         card = FALSE) {
  make_barcode_raw <- function() {
    #if multiple, then pull single "rep" image; consider pulling >1 and using patchwork, eg.
    if(length(input$content > 1)){
      gene_symbol <- sample(input$content, 1) #input$content[[1]]
    } else {
      gene_symbol <- input$content
    }

    #image type
    if(card == TRUE){image_type = "card"} else {image_type = "plot"}

    #check if exists
    url <- glue::glue("https://{Sys.getenv('AWS_BARCODE_BUCKET_ID')}.s3.amazonaws.com/{gene_symbol}_barcode_{image_type}.jpeg")
    status <- httr::GET(url) %>% httr::status_code()

    if(status == 200){
      return(url)
    } else {
      return(glue::glue("https://{Sys.getenv('AWS_BARCODE_BUCKET_ID')}.s3.amazonaws.com/error_barcode_{image_type}.jpeg"))
    }
  }
  #error handling
  tryCatch(make_barcode_raw(),
           error = function(e){
             message(e)
             if(card == TRUE){image_type = "card"} else {image_type = "plot"}
             glue::glue("https://{Sys.getenv('AWS_BARCODE_BUCKET_ID')}.s3.amazonaws.com/error_barcode_{image_type}.jpeg")
           })
}

## IDEOGRAM PLOT --------------------------------------------------------
#' Ideogram Plot
#'
#' Each point shows the location of the query gene(s) on human chromosomes.
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns an ideogram plot. If an error is thrown, then will return a bomb plot
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_ideogram(input = list(type = "gene", query = "ROCK1", content = "ROCK1"))
#' make_ideogram(input = list(type = "gene", query = "ROCK1", content = "ROCK1"), card = TRUE)
#' \dontrun{
#' make_ideogram(input = list(type = "gene", content = "ROCK1"))
#' }
make_ideogram <- function(data_gene_location = gene_location,
                          data_gene_chromosome = gene_chromosome,
                          input = list(),
                          card = FALSE) {
  make_ideogram_raw <- function() {
    #set baseline chromosome info
    chromosome_list <- dplyr::pull(data_gene_chromosome, id)
    chromosome_levels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
    chromosome_list <- forcats::fct_relevel(chromosome_list, chromosome_levels)

    #adjust by n so it bumps genes/bands off chromosome ends
    n <- 2000000 #chr1 is 250M

    #get some pq arms for the background images
    pq <-
      data_gene_chromosome %>%
      dplyr::mutate(centromere = (centromereposition_mbp * 1000000), #correct for n adjustment here
                    p = centromere + (n/2),
                    q = (basepairs+n) - p) %>% #correct for n adjustment here
      dplyr::select(chromosome_name = id, p, q) %>%
      dplyr::mutate(q_start = 0,
                    q_end = q,
                    p_start = q + 1,
                    p_end = q + p) %>%
      tidyr::pivot_longer(cols = c("p", "q"), names_to = "arm", values_to = "length") %>%
      dplyr::mutate(y_start = ifelse(arm == "q", q_start, p_start+n),
                    y_end = ifelse(arm == "q", q_end-n, p_end)) %>%
      dplyr::select(chromosome_name, arm, y_start, y_end)

    #get max lengths for left_join below to normalize transcript sites to maximum, so bands are plotted in the right direction
    max_lengths <-
      data_gene_chromosome %>%
      dplyr::select(chromosome_name = id, basepairs) %>%
      dplyr::mutate(basepairs = basepairs + n) #correct for n adjustment here

    #make some bands
    band_boundaries <-
      data_gene_location %>%
      dplyr::group_by(chromosome_name, band) %>%
      dplyr::summarize(abs_min = min(transcript_start),
                       abs_max = max(transcript_end),
                       length = abs_max - abs_min) %>%
      dplyr::ungroup() %>%
      dplyr::filter(chromosome_name %in% chromosome_list) %>%
      dplyr::left_join(max_lengths, by = "chromosome_name") %>%
      dplyr::mutate(min = basepairs - abs_min - (n/2), #chromosome positions count from the top, so need to reverse order for plotting
                    max = basepairs - abs_max - (n/2) #correct for n adjustment here
      ) %>%
      dplyr::mutate(alternate = dplyr::row_number(chromosome_name) %% 2)

    #adding if (is.null(gene_symbol)), gene_symbol=NULL allows me to generate a plot, even without gene_symbol data, #no data filters, #drop geom_point here
    #get location data for each gene query
    gene_loci <-
      data_gene_location %>%
      dplyr::filter(approved_symbol %in% input$content) %>%
      dplyr::select(approved_symbol, chromosome_name, transcript_start, transcript_end) %>%
      dplyr::left_join(max_lengths, by = "chromosome_name") %>% #need to flip the query gene start site around b/c chromosome positions count from the top, but ggplots from the bottom
      dplyr::mutate(start = basepairs - transcript_start - (n/2)) #correct for n adjustment here

    #error catching
    if(nrow(gene_loci) == 0){stop("gene not found")}

    #pull the choromosomes here, and then add a filter in the ggplot to plot them only
    chromosome_loci <-
      gene_loci %>%
      dplyr::distinct(chromosome_name) %>%
      dplyr::pull(.)

    #max lengths
    clip_height <-
      max_lengths %>%
      dplyr::filter(chromosome_name %in% chromosome_loci) %>%
      dplyr::slice_max(basepairs) %>%
      dplyr::pull(basepairs)

    ideogram_plot <-
      ggplot2::ggplot() +
      #background line fixes height
      ggplot2::geom_segment(data = pq %>% dplyr::filter(chromosome_name %in% chromosome_loci),
                            ggplot2::aes(x = chromosome_name, xend = chromosome_name, y = 0, yend = 260000000, alpha = 1),
                            color = "white", size = 1, lineend = "butt") +
      #background for black line
      ggplot2::geom_segment(data = pq %>% dplyr::filter(chromosome_name %in% chromosome_loci),
                            ggplot2::aes(x = chromosome_name, xend = chromosome_name, y = y_start, yend = y_end, alpha = 0.5),
                            color = "black", size = 7, lineend = "round") +
      #chromosome
      ggplot2::geom_segment(data = pq %>% dplyr::filter(chromosome_name %in% chromosome_loci),
                            ggplot2::aes(x = chromosome_name, xend = chromosome_name, y = y_start, yend = y_end),
                            color = "gray95", size = 6, lineend = "round") +
      #centromere
      ggplot2::geom_point(data = pq %>% dplyr::filter(chromosome_name %in% chromosome_loci,
                                                      arm == "q"),
                          ggplot2::aes(x = chromosome_name, y = y_end + n), #calculated 1/2 of distance between
                          color = "gray90", size = 5.5) +
      #pq bands
      ggplot2::geom_segment(data = band_boundaries %>% dplyr::filter(alternate == 1, chromosome_name %in% chromosome_loci),
                            ggplot2::aes(x = chromosome_name, xend = chromosome_name, y = min, yend = max),
                            size = 6.5, color = "black", alpha = 0.5) +
      #gene points + labels
      ggplot2::geom_point(data = gene_loci, ggplot2::aes(x = chromosome_name, y = start),  size = 6, color = ddh_pal_d(palette = "gene")(1), alpha = 1) +
      ggrepel::geom_text_repel(data = gene_loci, ggplot2::aes(x = chromosome_name, y = start, label = approved_symbol), nudge_x = .2, min.segment.length = 1, family = "Chivo") +
      ggplot2::scale_fill_manual(values = c("#FFFFFF", "#FFFFFF")) +
      ggplot2::labs(y = NULL) +
      theme_ddh(base_size = 16) +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = ggplot2::element_text(size = 12)) +
      #clip height
      ggplot2::coord_cartesian(ylim=c(0, clip_height)) +
      NULL

    if(card == TRUE){
      ideogram_plot <-
        ideogram_plot +
        ggplot2::labs(x = "") + #, title = "Gene Information", caption = "more ...") +
        NULL
    }
    return(ideogram_plot)
  }
  #error handling
  tryCatch(make_ideogram_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## SIZE PLOT --------------------------------------------------------
#' Protein Size Plot
#'
#' Mass compared to all protein masses. The colored strip visualizes the distribution of protein sizes. Each colored box is thus representing a decile of the full data. The triangle indicates where exactly the queried genes fall on this gradient of protein sizes.
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a proteinsize plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_proteinsize(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_proteinsize(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_proteinsize(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_proteinsize <- function(data_universal_proteins = universal_proteins,
                             input = list(),
                             card = FALSE) {
  make_proteinsize_raw <- function() {

    ### widths: distinct chunks or gradual strip?
    w <- seq(.1, 1, by = .1) ## 10 steps
    # w <- seq(.01, 1, by = .01) ## 100 steps ~ almost gradient
    # w <- c(0.5, 0.95, 1) ## threshold values: 50%, 95%, 100% (or smt similar)
    colors <- ddh_pal_c(palette = "protein")(length(w))

    ## sort alphabetically
    #gene_symbol <- sort(gene_symbol)

    ## sort by mass
    gene_symbol <-
      data_universal_proteins %>%
      dplyr::filter(gene_name %in% input$content) %>%
      dplyr::arrange(mass) %>%
      dplyr::pull(gene_name)

    moreThanThree <- length(gene_symbol) > 3
    moreThanFour<- length(gene_symbol) > 4

    last_gene <- gene_symbol[length(gene_symbol)]
    last_mass <-
      data_universal_proteins %>%
      dplyr::filter(gene_name == last_gene) %>%
      dplyr::pull(mass)

    make_mass_strip <- function(data = data_universal_proteins,
                                var,
                                card_var = card,
                                max_mass = last_mass) {

      selected <-
        data %>%
        dplyr::filter(gene_name %in% var)

      if(card_var == TRUE){ #this solves for a single case
        triangle_size <- 10
        triangle_y <- 1.025
        mass_y <- 1.04
        gene_y <- 1.05
      } else {
        triangle_size <- 6
        triangle_y <- 1.06
        mass_y <- 1.085
        gene_y <- 1.12
      }

      base_plot <-
        data %>%
        ggplot2::ggplot(ggplot2::aes(x = mass)) +
        ## draw distribution as colored strip
        ggdist::stat_interval(
          ggplot2::aes(y = 1),
          .width = w, size = 10
        ) +
        ## line locator
        ggplot2::geom_linerange(
          data = selected, ggplot2::aes(ymin = .9, ymax = 1.03),
          color = "white", size = .8
        ) +
        ## add triangular locator symbol
        ggplot2::geom_point(
          data = selected, ggplot2::aes(y = triangle_y),
          shape = 6, size = triangle_size, stroke = 1
        ) +
        ## add gene name label
        {if(!card_var)ggplot2::geom_text( #| !moreThanFour
          data = selected, ggplot2::aes(y = gene_y, label = gene_name),
          family = "Chivo", size = 5.2, vjust = 0
        )} +
        ## add mass label
        {if(!card_var)ggplot2::geom_text( # | !moreThanThree
          data = selected, ggplot2::aes(y = mass_y, label = paste(mass, "kDa")),
          family = "Chivo", size = 4.3, vjust = 0
        )} +
        ggplot2::scale_x_continuous(limits = c(0, max_mass*2),
                                    labels = function(x) paste(x, "kDa")) +
        ggplot2::coord_cartesian(ylim = c(.95, 1.2)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::scale_color_manual(values = colors, guide = "none") +
        theme_ddh() +
        ggplot2::theme_void() +
        NULL


      if (var == last_gene) {
        base_plot <-
          base_plot +
          ggplot2::labs(x = "Protein Size") +
          ggplot2::theme(axis.text.x = ggplot2::element_text(family = "Roboto Slab", size = 14),
                         axis.title.x = ggplot2::element_text(family = "Nunito Sans", size = 18,
                                                              margin = ggplot2::margin(t = 12, b = 12)))
      }

      return(base_plot)
    }

    glist <- lapply(gene_symbol, make_mass_strip, data = data_universal_proteins)
    plot_complete <- patchwork::wrap_plots(glist, nrow = length(glist))

    if(card == TRUE){
      plot_complete <- make_mass_strip(var = gene_symbol)

      mass <- #get mass to center clipping in next step
        data_universal_proteins %>%
        dplyr::filter(gene_name %in% input$content) %>%
        dplyr::pull(mass) %>%
        median(na.rm = TRUE)

      plot_complete <-
        plot_complete +
        ggplot2::coord_cartesian(ylim = c(.9, 1.1), xlim = c(mass - mass*.3, mass + mass*.3)) +
        # plot_annotation(
        #   title = "Size Information",
        #   caption = "more ...") +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank()
        )
    }

    return(plot_complete)
  }
  #error handling
  tryCatch(make_proteinsize_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## SEQUENCE PLOT --------------------------------------------------------------------
#' Sequence Plot
#'
#' \code{make_sequence} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a sequence plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a sequence plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_sequence(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_sequence(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_sequence(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_sequence <- function(data_universal_proteins = universal_proteins,
                          input = list(),
                          card = FALSE) {
  make_sequence_raw <- function() {
    sequence_string <-
      data_universal_proteins %>%
      dplyr::filter(gene_name %in% input$content) %>%
      dplyr::pull("sequence") %>%
      stringr::str_sub(start = 1L, end = 100L) %>%
      stringr::str_pad(100, "right")

    #compose sequence vec
    starting_numbers <- seq(from = 1, to = 100, by=10)
    sequence_num <- vctrs::vec_rep_each(starting_numbers, length(sequence_string))
    #sequence_sub_vec <- vctrs::vec_rep(1:length(sequence_string), 10) #put 10 here for the case of one sequence
    sequence_vec <- vctrs::vec_rep(sequence_string[1:length(sequence_string)], 10)

    get_subsequence <- function(seq_string, start_seq){
      subsequence <- stringr::str_sub(seq_string,
                                      start = start_seq,
                                      end = start_seq + 9)
      return(subsequence)
    }
    #get_subsequence(seq_string = sequence_string, seq_sub = 1, start_seq = 1)

    sequence_complete <-
      purrr::map2(.x = sequence_vec,
                  .y = sequence_num,
                  .f = get_subsequence)

    sequence_font <- "Roboto Slab"
    sequence_size <- 10

    plot_complete <-
      ggplot2::ggplot() +
      ## annotate with text
      ggplot2::annotate("text", x = 0, y = 10, label = sequence_complete[1], family = sequence_font, size = sequence_size) +
      ggplot2::annotate("text", x = 0, y = 9, label = sequence_complete[2], family = sequence_font, size = sequence_size) +
      ggplot2::annotate("text", x = 0, y = 8, label = sequence_complete[3], family = sequence_font, size = sequence_size) +
      ggplot2::annotate("text", x = 0, y = 7, label = sequence_complete[4], family = sequence_font, size = sequence_size) +
      ggplot2::annotate("text", x = 0, y = 6, label = sequence_complete[5], family = sequence_font, size = sequence_size) +
      ggplot2::annotate("text", x = 0, y = 5, label = sequence_complete[6], family = sequence_font, size = sequence_size) +
      ggplot2::annotate("text", x = 0, y = 4, label = sequence_complete[7], family = sequence_font, size = sequence_size) +
      ggplot2::annotate("text", x = 0, y = 3, label = sequence_complete[8], family = sequence_font, size = sequence_size) +
      ggplot2::annotate("text", x = 0, y = 2, label = sequence_complete[9], family = sequence_font, size = sequence_size) +
      ggplot2::annotate("text", x = 0, y = 1, label = sequence_complete[10], family = sequence_font, size = sequence_size) +
      ggplot2::coord_cartesian(ylim = c(0.5, 10.5), clip = 'on') +
      ggplot2::theme_void() +
      NULL

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        ggplot2::labs(x = "") +#, title = "Sequence Information", caption = "more ...") +
        NULL
    }

    return(plot_complete)
  }
  #error handling
  tryCatch(make_sequence_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## PROTEIN DOMAIN PLOT --------------------------------------------------------
#' Protein Domain Plot
#'
#' Rectangles represent the locations and size of named protein domains, while black shaped elements represent PTMs. Horizontal line(s) indicate the length of one or more selected proteins.
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a protein domain plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_protein_domain(input = list(content = "ROCK2"), dom_var = "Protein kinase", ptm_var = "N-acetylserine")
#' make_protein_domain(input = list(content = c("ROCK1", "ROCK2")), dom_var = "Protein kinase", ptm_var = "N-acetylserine")
#' \dontrun{
#' make_protein_domain(input = list(content = "ROCK2"), dom_var = "Protein kinase", ptm_var = "N-acetylserine")
#' }
make_protein_domain <- function(input = list(),
                                data_gene_protein_domains = gene_protein_domains,
                                dom_var = NULL,
                                ptm_var = NULL) {

  make_protein_domain_raw <- function() {

    gene_symbol <-
      data_gene_protein_domains %>%
      dplyr::filter(gene_name %in% input$content) %>%
      dplyr::pull(gene_name) %>%
      unique()

    lengths_data <-
      data_gene_protein_domains %>%
      dplyr::filter(gene_name %in% input$content) %>%
      dplyr::filter(!duplicated(gene_name))

    plot_data <-
      data_gene_protein_domains %>%
      dplyr::filter(gene_name %in% input$content) %>%
      dplyr::mutate(order = 1)

    prots_dr <-
      plot_data %>%
      dplyr::filter(type %in% c("DOMAIN", "REGION")) %>%
      tidyr::drop_na(description)

    prots_ptm <-
      plot_data %>%
      dplyr::filter(category == "PTM") %>%
      tidyr::drop_na(description)

    # Plot
    base_plot <-
      ggplot2::ggplot() +
      ggplot2::geom_segment(data = plot_data,
                            ggplot2::aes(x = 1,
                                         xend = seq_len,
                                         y = order,
                                         yend = order),
                            colour = "black",
                            size = 1.5,
                            lineend = "round")
    #DOMAINS
    if(nrow(prots_dr) != 0) {
      base_plot <-
        base_plot +
        ggplot2::geom_rect(data = prots_dr %>%
                             dplyr::filter(description %in% dom_var),
                           ggplot2::aes(xmin = begin,
                                        xmax = end,
                                        ymin = order - 1,
                                        ymax = order + 1,
                                        fill = description),
                           alpha = 0.9,
                           color = "black")
    }

    #PTMS
    if(nrow(prots_ptm) != 0) {
      base_plot <-
        base_plot +
        ggplot2::geom_point(data = prots_ptm %>%
                              dplyr::filter(description %in% ptm_var),
                            ggplot2::aes(x = begin,
                                         y = order + 3,
                                         shape = description),
                            size = 8,
                            alpha = 1,
                            color = "black") +
        ggplot2::geom_segment(data = prots_ptm %>%
                                dplyr::filter(description %in% ptm_var),
                              ggplot2::aes(x = begin,
                                           xend = begin,
                                           y = order,
                                           yend = order + 3),
                              size = 0.5,
                              linetype = 2,
                              color = "black") +
        ggplot2::scale_y_continuous(limits=c(0,4.1)) + #should be a touch larger than yend
        NULL
    }

    wrapped_labels <- function(.seq_len) {
      fun <- function(labels) {
        labels <- ggplot2::label_value(labels, multi_line = TRUE)
        lapply(labels, function(x) {
          x <- paste0(x, " (", .seq_len, " amino acids)")
          vapply(x, paste, character(1), collapse = "\n")
        })
      }
      structure(fun, class = "labeller")
    }

    base_plot <-
      base_plot +
      ggplot2::labs(x = NULL,
                    y = NULL) +
      ggplot2::facet_wrap(~ gene_name, ncol = 1,
                          labeller = wrapped_labels(.seq_len = lengths_data$seq_len)
      ) +
      theme_ddh(base_size = 16) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(family = "Chivo", hjust = 0.5, size = 15),
        plot.title = ggplot2::element_text(family = "Chivo", hjust = 0.5, size = 15),
        text = ggplot2::element_text(family = "Nunito Sans"),
        legend.text = ggplot2::element_text(size = 15),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
      ) +
      NULL

    plot_complete <-
      base_plot +
      ggplot2::theme(axis.text.x = ggplot2::element_text(family = "Roboto Slab", size = 14),
                     axis.title.x = ggplot2::element_text(family = "Nunito Sans", size = 18,
                                                          margin = ggplot2::margin(t = 12, b = 12)),
                     legend.position = "right") +
      scale_fill_ddh_d(palette = "protein")

    return(plot_complete)

  }

  #error handling
  tryCatch(make_protein_domain_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## AA RADIAL PLOT -------------------------------------------------------------
#' Amino Acid Radial Plot
#'
#' Amino acid signature/s (percentage of each amino acid in a protein) of the queried gene/clusters versus the mean amino acid signature of all the other proteins in the dataset (N = 20375).
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a radial plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#' @export
#' @examples
#' make_radial(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_radial(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), cluster = TRUE)
#' make_radial(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_radial(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_radial <- function(data_gene_signatures = gene_signatures,
                        data_gene_signature_clusters = gene_signature_clusters,
                        input = list(),
                        relative = TRUE,
                        cluster = FALSE,
                        barplot = FALSE,
                        card = FALSE) {
  make_radial_raw <- function() {
    if(cluster) {
      query_clust <-
        data_gene_signature_clusters %>%
        dplyr::filter(gene_name %in% input$content) %>%
        dplyr::pull(clust) %>%
        unique()

      signature_cluster_means_prep <-
        data_gene_signatures %>%
        dplyr::select(uniprot_id, A:Y) %>%
        dplyr::inner_join(data_gene_signature_clusters, by = "uniprot_id") %>%
        dplyr::select(-gene_name, -X1, -X2, -member_prob)

      signature_cluster_means_query <-
        signature_cluster_means_prep %>%
        dplyr::mutate(clust = as.numeric(as.character(clust)),
                      clust = as.factor(ifelse(clust %in% as.numeric(as.character(query_clust)), clust, "Mean"))) %>%
        dplyr::group_by(clust) %>%
        dplyr::summarise_if(is.numeric, list(mean = mean)) %>%
        tidyr::pivot_longer(cols = -clust) %>%
        dplyr::mutate(name = stringr::str_remove(name, "_mean"),
                      clust = as.factor(ifelse(clust == "Mean", "Mean",
                                               paste0("Cluster ", as.numeric(as.character(clust)))))
        )

      if(relative){
        signature_cluster_means_query <-
          signature_cluster_means_query %>%
          tidyr::pivot_wider(id_cols = clust, values_fn = mean) %>%
          dplyr::arrange(factor(clust, levels = "Mean")) %>%
          dplyr::mutate(across(A:Y, ~ . / .[1])) %>%
          tidyr::pivot_longer(cols = -clust)
      }

      if(length(unique(signature_cluster_means_query$clust)) == 1) {
        stop("Unable to cluster this protein by its amino acid sequence.")
      }

    } else {
      signature_cluster_means_prep <-
        data_gene_signatures %>%
        dplyr::select(gene_name, A:Y)

      signature_cluster_means_query <-
        signature_cluster_means_prep %>%
        dplyr::filter(!gene_name %in% input$content) %>%
        dplyr::rename(Mean = gene_name) %>%
        dplyr::summarise_if(is.numeric, list(mean = mean)) %>%
        tibble::rownames_to_column("clust") %>%
        tidyr::pivot_longer(cols = -clust) %>%
        dplyr::mutate(name = stringr::str_remove(name, "_mean"),
                      clust = "Mean")

      signature_gene_query <-
        signature_cluster_means_prep %>%
        dplyr::filter(gene_name %in% input$content) %>%
        tidyr::pivot_longer(cols = -gene_name) %>%
        dplyr::rename(clust = gene_name)

      signature_cluster_means_query <-
        dplyr::bind_rows(signature_cluster_means_query,
                         signature_gene_query)

      if(relative){
        signature_cluster_means_query <-
          signature_cluster_means_query %>%
          tidyr::pivot_wider(id_cols = clust, values_fn = mean) %>%
          dplyr::arrange(factor(clust, levels = "Mean")) %>%
          dplyr::mutate(across(A:Y, ~ . / .[1])) %>%
          tidyr::pivot_longer(cols = -clust)
      }
    }

    # set colors -1 to assign specific color to "Mean"
    colors_raw <- ddh_pal_d(palette = "protein")(length(unique(signature_cluster_means_query$clust))-1)
    names(colors_raw) <- unique(signature_cluster_means_query$clust)[unique(signature_cluster_means_query$clust) != "Mean"]
    mean_color <- "gray48"
    names(mean_color) <- "Mean"
    colors_radial <- c(colors_raw, mean_color)

    # barplot fill
    mean_fill_bar <- "gray94"
    names(mean_fill_bar) <- "Mean"
    colors_bar <- c(colors_raw, mean_fill_bar)

    #relative label
    if(relative == TRUE){
      y_label = "Relative AA Frequency"
    } else {
      y_label = "AA Frequency (%)"
    }
    # RADIAL/BAR PLOT
    plot_complete <-
      ggplot2::ggplot(signature_cluster_means_query,
                                     ggplot2::aes(x = forcats::fct_inorder(name),
                                                  y = value,
                                                  group = clust,
                                                  color = clust
                                     )) +
      {if(!barplot)ggplot2::geom_point(alpha = 0.8, show.legend = FALSE)} +
      {if(!barplot)ggplot2::geom_polygon(fill = NA)} +
      {if(barplot & relative)ggplot2::geom_col(data = signature_cluster_means_query %>%
                                                 dplyr::filter(clust != "Mean"),
                                               ggplot2::aes(x = reorder(name, -value),
                                                            y = value,
                                                            group = clust,
                                                            color = clust,
                                                            fill = clust),
                                               position = "dodge2")} +
      {if(barplot & !relative)ggplot2::geom_col(aes(fill = clust), position = "dodge2")} +
      {if(barplot & relative)ggplot2::geom_hline(yintercept = 1, color = "gray48")} +
      ggplot2::labs(y = y_label,
                    x = ggplot2::element_blank()) +
      {if(!barplot)ggplot2::coord_polar()} +
      ggplot2::scale_color_manual(values = colors_radial) +
      {if(barplot)ggplot2::scale_fill_manual(values = colors_bar)} +
      ## theme changes
      theme_ddh(base_size = 16) +
      {if(!barplot)ggplot2::theme_minimal()} +
      {if(!barplot)ggplot2::theme(
        text = ggplot2::element_text(family = "Nunito Sans"),
        legend.position = "top",
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 15),
        axis.line.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(family = "Roboto Slab", size = 16),
        axis.text.y = ggplot2::element_text(size = 16, color = "grey30"),
        axis.title = ggplot2::element_text(size = 16)
      )} +
      {if(barplot)ggplot2::theme(
        text = ggplot2::element_text(family = "Nunito Sans"),
        legend.position = "top",
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 15),
        axis.text = ggplot2::element_text(family = "Roboto Slab"),
        axis.ticks.x = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank()
      )} +
      NULL

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        ggplot2::labs(x = "", y = "") + #, title = "Signature Information", caption = "more ...") +
        ggplot2::theme(
          text = ggplot2::element_text(family = "Nunito Sans"),
          legend.position = "none",
          legend.title = ggplot2::element_blank(),
          legend.text = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(family = "Roboto Slab"),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_blank()
        )
    }

    return(plot_complete)
  }

  #error handling
  tryCatch(make_radial_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## AA BAR PLOT -------------------------------------------------------------
#' Amino Acid Bar Plot
#'
#' Amino acid signature/s (percentage of each amino acid in a protein) of the queried gene/clusters versus the mean amino acid signature of all the other proteins in the dataset (N = 20375).
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a bar plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#' @export
#' @examples
#' make_radial_bar(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_radial_bar(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), cluster = TRUE)
#' make_radial_bar(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_radial_bar(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_radial_bar <- function(input = list(),
                            relative = TRUE,
                            cluster = FALSE,
                            card = FALSE) {
  make_radial_bar_raw <- function() {
    plot_complete <- make_radial(input = input,
                                 relative = relative,
                                 cluster = cluster,
                                 barplot = TRUE,
                                 card = card)

    return(plot_complete)
  }

  #error handling
  tryCatch(make_radial_bar_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## UMAP PLOT --------------------------------------------------------
#' UMAP Plot
#'
#' Amino acid signature (percentage of each amino acid in a protein) UMAP embeddings (2D) colored by the cluster to which they belong (N = 20375).
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a umap plot plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_umap_plot(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_umap_plot(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_umap_plot(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_umap_plot <- function(data_gene_signature_clusters = gene_signature_clusters,
                           input = list(),
                           show_subset = FALSE,
                           labels = FALSE) {
  make_umap_plot_raw <- function() {

    data_proteins_clean <-
      data_gene_signature_clusters %>%
      # dplyr::filter(clust != 0) %>%
      dplyr::mutate(clust = paste0("Cluster ", clust))

    query_clust <-
      data_proteins_clean %>%
      dplyr::filter(gene_name %in% input$content) %>%
      dplyr::pull(clust) %>%
      unique()

    cluster_genes <-
      data_proteins_clean %>%
      dplyr::filter(clust %in% query_clust)

    colors <- ddh_pal_c(palette = "protein")(length(query_clust))

    # UMAP PLOT
    plot_complete <- ggplot2::ggplot() +
      {if(!show_subset)ggplot2::geom_point(data = data_proteins_clean %>%
                                             dplyr::filter(!uniprot_id %in% cluster_genes$uniprot_id),
                                           ggplot2::aes(X1, X2), size = 0.8, color = "grey80")} +
      ggplot2::geom_point(data = data_proteins_clean %>%
                            dplyr::filter(uniprot_id %in% cluster_genes$uniprot_id),
                          ggplot2::aes(X1, X2, color = clust), size = 0.8) +
      {if(labels)ggrepel::geom_label_repel(data = cluster_genes %>%
                                             dplyr::filter(gene_name %in% input$content),
                                           ggplot2::aes(X1, X2, label = gene_name))} +
      ggplot2::labs(x = "UMAP 1",
                    y = "UMAP 2") +
      ggplot2::scale_color_manual(
        values = rep(colors, length.out =
                       nrow(data_proteins_clean %>%
                              dplyr::filter(uniprot_id %in% cluster_genes$uniprot_id)))
      ) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3))) +
      ## theme changes
      theme_ddh() +
      ggplot2::theme(
        text = ggplot2::element_text(family = "Nunito Sans"),
        legend.position = "top",
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 15),
        axis.text = ggplot2::element_text(family = "Roboto Slab", size = 16),
        axis.text.y = ggplot2::element_text(size = 16, color = "grey30"),
        axis.title = ggplot2::element_text(size = 16)
      ) +
      NULL

    return(plot_complete)
  }

  #error handling
  tryCatch(make_umap_plot_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## CLUSTER ENRICHMENT PLOT --------------------------------------------------------
#' Cluster Enrichment Plot
#'
#' Enriched GO terms (BP, MF, and CC) by all genes in the selected amino acid signature cluster. The x-axis shows the number of genes in the cluster that belong to each term while the color scale represents the p-values of each enriched term.
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a cluster enrich plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_cluster_enrich(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_cluster_enrich(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_cluster_enrich(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_cluster_enrich <- function(data_gene_signature_clusters = gene_signature_clusters,
                                data_gene_signature_cluster_enrichment = gene_signature_cluster_enrichment,
                                input = list(),
                                ontology = "BP",
                                num_terms = 20) {

  make_cluster_enrich_plot_raw <- function() {

    plot_data <- make_clustering_enrichment_table(data_gene_signature_clusters,
                                                  data_gene_signature_cluster_enrichment,
                                                  input = input,
                                                  ontology = ontology)

    plot_complete <-
      plot_data %>%
      dplyr::arrange(pvalue) %>%
      dplyr::slice(1:num_terms) %>%
      ggplot2::ggplot(ggplot2::aes(x = Count,
                                   y = reorder(substr(paste0(Description, " (", ID, ")"), 1, 40), Count),
                                   fill = pvalue)) +
      ggplot2::geom_col() +
      ggplot2::labs(x = "Gene Count",
                    y = NULL) +
      scale_fill_ddh_c(palette = "protein", reverse = TRUE) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(barheight = ggplot2::unit(6, "lines"),
                                                     barwidth = ggplot2::unit(.6, "lines"),
                                                     reverse = TRUE)) +
      theme_ddh() +
      ggplot2::theme(
        title = ggplot2::element_blank(),
        text = ggplot2::element_text(family = "Nunito Sans"),
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 16),
        legend.text = ggplot2::element_text(size = 16),
        axis.text = ggplot2::element_text(family = "Roboto Slab"),
        axis.ticks.x = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank()
      ) +
      NULL

    return(plot_complete)
  }
  #error handling
  tryCatch(make_cluster_enrich_plot_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## STRUCTURE PLOT --------------------------------------------------------------------
#' Protein Structure Plot
#'
#' Alpha Fold predicted structure rendered in a ribbon diagram.
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a structure plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_structure(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_structure(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_structure(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_structure <- function(input = list(),
                           card = FALSE){
  make_structure_raw <- function(){
    #if multiple, then pull single "rep" image; consider pulling >1 and using patchwork, eg.
    if(length(input$content > 1)){
      gene_symbol <- sample(input$content, 1) #input$content[[1]]
    } else {
      gene_symbol <- input$content
    }

    #image type
    if(card == TRUE){image_type = "card"} else {image_type = "plot"}

    #check if exists
    url <- glue::glue("https://{Sys.getenv('AWS_PROTEINS_BUCKET_ID')}.s3.amazonaws.com/{gene_symbol}_structure_{image_type}.jpg")
    status <- httr::GET(url) %>% httr::status_code()

    if(status == 200){
      return(url)
    } else {
      return(glue::glue("https://{Sys.getenv('AWS_PROTEINS_BUCKET_ID')}.s3.amazonaws.com/error_structure_{image_type}.jpg"))
    }
  }
  #error handling
  tryCatch(make_structure_raw(),
           error = function(e){
             message(e)
             if(card == TRUE){image_type = "card"} else {image_type = "plot"}
             glue::glue("https://{Sys.getenv('AWS_PROTEINS_BUCKET_ID')}.s3.amazonaws.com/error_structure_{image_type}.jpg")
           })
}

## 3D STRUCTURE PLOT --------------------------------------------------------------------
#' Protein 3D Structure Plot
#'
#' Protein 3D predicted structure rendered in an interactive ribbon diagram.
#'
#' @param input Expecting a list containing type and content variable.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_structure3d(input = list(type = 'gene', content = 'ROCK1'))
#' make_structure3d(input = list(type = 'gene', content = 'ROCK2'))
#' \dontrun{
#' make_structure3d(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_structure3d <- function(data_gene_uniprot_pdb_table = gene_uniprot_pdb_table,
                             data_universal_proteins = universal_proteins,
                             gene_symbol = NULL,
                             pdb_id = NULL,
                             input = list(),
                             color = FALSE,
                             ribbon = FALSE,
                             selection = FALSE,
                             resi = 1:10,
                             chain = "A",
                             resn = NULL,
                             invert = NULL,
                             elem = NULL
) {
  make_structure3d_raw <- function() {
    #because this fun doesn't take multi-gene queries
    if(is.null(gene_symbol)) {
      if(length(input$content) == 1) {
        gene_symbol <- input$content
      } else {
        gene_symbol <- input$content[1]
      }
    }
    #check  to see if pdb file exists
    pdb_path <- load_pdb(input = input)

    if(!is.null(pdb_path) & is.null(pdb_id)) {
      plot_data <- bio3d::read.pdb(pdb_path)
    } else {
      plot_data <-
        data_universal_proteins %>%
        dplyr::filter(gene_name %in% gene_symbol) %>%
        dplyr::left_join(data_gene_uniprot_pdb_table, by = c("uniprot_id" = "uniprot")) %>%
        tidyr::unnest(data) %>%
        {if (is.null(pdb_id)) dplyr::slice(., 1) else dplyr::filter(., pdb == pdb_id)} %>%
        dplyr::pull(pdb) %>%
        r3dmol::m_fetch_pdb()
    }

    plot_complete <-
      r3dmol::r3dmol() %>%
      r3dmol::m_add_model(data = plot_data, format = "pdb") %>%
      r3dmol::m_center()

    if(color) {
      plot_complete <- plot_complete %>%
        r3dmol::m_set_style(style = r3dmol::m_style_cartoon(color = "spectrum", ribbon = ribbon)) %>%
        r3dmol::m_center()
    } else {
      plot_complete <- plot_complete %>%
        r3dmol::m_set_style(style = r3dmol::m_style_cartoon(ribbon = ribbon)) %>%
        r3dmol::m_center()
    }

    if(selection) {
      plot_complete <- plot_complete %>%
        r3dmol::m_add_style(
          style = c(
            r3dmol::m_style_stick(),
            r3dmol::m_style_sphere(scale = 0.3)
          ),
          sel = r3dmol::m_sel(resi =  eval(parse(text = resi)),
                              chain = chain,
                              resn = resn,
                              invert = invert,
                              elem = elem)
        ) %>%
        m_zoom_to(sel = m_sel(resi = eval(parse(text = resi)),
                              chain = chain)) %>%
        r3dmol::m_center()
    }

    return(plot_complete)
  }

  #error handling
  tryCatch(make_structure3d_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## PUBMED PLOT ------------------------------------------
#' Pubmed Plot
#'
#' \code{make_pubmed} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a pubmed plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a pubmed plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_pubmed(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_pubmed(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' make_pubmed(input = list(type = 'compound', content = 'aspirin'))
#' \dontrun{
#' make_pubmed(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_pubmed <- function(data_universal_pubmed = universal_pubmed,
                        input = list(),
                        card = FALSE) {
  make_pubmed_raw <- function() {
    plot_data <-
      data_universal_pubmed %>%
      dplyr::filter(name %in% input$content) %>%
      tidyr::unnest(data) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(name, year) %>%
      dplyr::summarize(n = dplyr::n()) %>%
      dplyr::mutate(cumsum = cumsum(n)) %>%
      dplyr::arrange(year)

    plot_max <-
      plot_data %>%
      dplyr::filter(!is.na(year)) %>%
      dplyr::filter(cumsum == max(cumsum))

    plot_step <-
      ggplot2::ggplot(plot_data) +
      ggplot2::geom_step(
        ggplot2::aes(x = year,
                     y = cumsum,
                     group = name,
                     color = forcats::fct_reorder2(name, year, cumsum)),
        size = 1.2
      ) +
      ggplot2::coord_cartesian(clip = "off") + #allows points & labels to fall off plotting area
      ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(4),
                                  expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(.01, .01), limits = c(0, max(plot_max$cumsum))) +
      scale_color_ddh_d(palette = input$type) +
      ggplot2::labs(x = "Year of Publication", y = "Cumulative Sum") +
      theme_ddh(base_size = 16,
                margin = 20) +
      ggplot2::theme(
        #text = element_text(family = "Nunito Sans"),
        #axis.text.y = element_text(family = "Roboto Slab"),
        panel.grid.major.y = ggplot2::element_line(color = "grey75", size = .5, linetype = "15")
      ) +
      NULL

    ## only add labels + white space if several genes queried
    if (length(unique(plot_data$name)) > 1) {
      plot_complete <-
        plot_step +
        {if(card == FALSE)
          ggrepel::geom_text_repel(
            data = plot_max,
            ggplot2::aes(x = year,
                         y = cumsum,
                         #label = name),
                         label = paste0(name, " (", cumsum, ")")),
            size = 5.5,
            hjust = 0,
            direction = "y",
            nudge_x = 1, # bump 2 years to the right
            xlim = c(max(plot_max$year) + 2, Inf), # no constraints on right side for labels
            max.overlaps = Inf, #5
            family = "Roboto Slab",
            segment.color = "grey65",
            segment.size = 0.9,
            min.segment.length = 0, # always draw segment
            force = .5,
            segment.curvature = -0.15,
            segment.ncp = 3,
            segment.angle = 90,
            segment.inflect = FALSE,
            box.padding = .2
          )
        } +
        ggplot2::geom_point(
          data = plot_max,
          ggplot2::aes(x = year,
                       y = cumsum,
                       group = name,
                       color = forcats::fct_reorder2(name, year, cumsum)),
          size = 4, shape = 21, fill = "white", stroke = 2
        ) +
        ggplot2::theme(
          plot.margin = ggplot2::margin(7, 180, 7, 7) #adds margin to right side of graph for label #adds margin to right side of graph for label
        )
    } else {
      plot_complete <-
        plot_step +
        ggplot2::geom_point(
          data = plot_max,
          ggplot2::aes(x = year,
                       y = cumsum,
                       group = name,
                       color = forcats::fct_reorder2(name, year, cumsum)),
          size = 4, shape = 21, fill = "white", stroke = 2
        )
    }

    if(length(input$content) == 1){ #Fix this when data()$content is updated (need to be generic to handle multiple data types)
      plot_complete  <-
        plot_complete +
        ggplot2::labs(y = "Cumulative Publications",
                      color = "") +
        ggplot2::guides(color = "none")
    } else {
      plot_complete <-
        plot_complete +
        ggplot2::labs(y = "Cumulative Publications",
                      color = "Query")
    }

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        gplot2::labs(x = "") +
        ggplot2::theme(plot.margin = ggplot2::margin(5, 10, 5, 5),
                       legend.position="none",
                       axis.text.x=element_blank()) +
        NULL
    }

    return(plot_complete)
  }
  #error handling
  tryCatch(make_pubmed_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## CELL ANATOGRAM -------------------------------------------------------------------
#' Cellanatogram Plot
#'
#' \code{make_cellanatogram} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a cellanatogram plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a cellanatogram plot. If an error is thrown, then will return a bomb plot.
#'
#' @import gganatogram
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_cellanatogram(input = list(type = "gene", content = c("ROCK2")))
#' make_cellanatogram(input = list(type = "gene", content = c("ROCK2")), card = TRUE)
#' make_cellanatogram(input = list(type = "gene", query = "ROCK1", content = c("ROCK1", "ROCK2")))
#' \dontrun{
#' make_cellanatogram(input = list(type = 'gene', content = 'ROCK2'))
#' }
make_cellanatogram <- function(data_gene_subcell = gene_subcell,
                               input = list(),
                               card = FALSE) {
  make_cellanatogram_raw <- function() {
    plot_data <-
      data_gene_subcell %>%
      dplyr::filter(gene_name %in% input$content) %>%
      dplyr::select(organ, type, colour, value) %>%
      tidyr::drop_na() %>%
      dplyr::group_by(organ) %>%
      dplyr::summarise(type = type[1], value = mean(value)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(value))

    plot_complete <-
      plot_data %>%
      gganatogram::gganatogram(outline = TRUE,
                               fillOutline = "grey95",
                               organism = "cell",
                               fill = "value") +
      ggplot2::theme_void(base_size = 14) +
      ggplot2::theme(
        text = ggplot2::element_text(family = "Nunito Sans"),
        plot.margin = ggplot2::margin(5, 10, 5, 5)
      ) +
      ggplot2::coord_fixed() +
      scale_fill_ddh_c(palette = "protein") +
      ggplot2::labs(fill = "Expression") +
      NULL

    if(length(input$content) == 1){
      plot_complete  <-
        plot_complete +
        ggplot2::guides(fill = "none")
    }

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        ggplot2::labs(x = "") +
        ggplot2::guides(fill = "none") +
        NULL
    }
    return(plot_complete)
  }
  #error handling
  tryCatch(make_cellanatogram_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

#' Cell Anatogram Facet
#'
#' @import gganatogram
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_cellanatogramfacet(input = list(type = "gene", content = c("BRCA1", "ROCK2")))
make_cellanatogramfacet <- function(data_gene_subcell = gene_subcell,
                                    input = list()) {
  make_cellanatogramfacet_raw <- function() {
    plot_data <-
      data_gene_subcell %>%
      dplyr::filter(gene_name %in% input$content) %>%
      dplyr::select(gene_name, organ, type, colour, value) %>%
      tidyr::drop_na() %>%
      dplyr::group_by(organ) %>%
      dplyr::summarise(type = gene_name, value = mean(value)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(value))

    plot_complete <-
      plot_data %>%
      gganatogram::gganatogram(outline = TRUE,
                               fillOutline = "grey95",
                               organism = "cell",
                               fill = "value") +
      ggplot2::theme_void(base_size = 14) +
      ggplot2::theme(
        text = ggplot2::element_text(family = "Nunito Sans"),
        plot.margin = ggplot2::margin(5, 10, 5, 5)
      ) +
      ggplot2::coord_fixed() +
      scale_fill_ddh_c(palette = "protein") +
      ggplot2::labs(fill = "Expression") +
      ggplot2::facet_wrap(~ type, ncol = 3, drop = FALSE) +
      ggplot2::theme(panel.spacing.y = ggplot2::unit(1.2, "lines"),
                     strip.text = ggplot2::element_text(size = 18, family = "Roboto Slab", face = "plain")) +
      NULL

    return(plot_complete)
  }
  #error handling
  tryCatch(make_cellanatogramfacet_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## HUMAN BODY ANATOGRAMS --------------------------------------------------------
#' Female Anatogram Plot
#'
#' \code{make_female_anatogram} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a female anatogram plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a female anatogram plot. If an error is thrown, then will return a bomb plot.
#'
#' @import gganatogram
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_female_anatogram(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_female_anatogram(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_female_anatogram(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_female_anatogram <- function(data_gene_female_tissue = gene_female_tissue,
                                  data_gene_male_tissue = gene_male_tissue,
                                  anatogram = "female",
                                  input = list(),
                                  card = FALSE) {
  make_female_anatogram_raw <- function() {
    if(anatogram == "female"){
      data_tissue <- data_gene_female_tissue
    } else if (anatogram == "male") {
      data_tissue <- data_gene_male_tissue
    } else {
      print("Declare your anatogram type")
    }

    body_data <-
      data_tissue %>%
      dplyr::filter_all(dplyr::any_vars(gene_name %in% input$content)) %>%
      dplyr::filter(!is.na(type)) %>%
      dplyr::arrange(dplyr::desc(-value))

    break_points <- unname(quantile(body_data$value))

    if(length(input$content) > 1){
      body_data <-
        body_data %>%
        dplyr::group_by(organ) %>%
        dplyr::mutate(value=sum(value)) %>%
        dplyr::arrange(dplyr::desc(-value))
    }

    body_plot <-
      body_data %>%
      gganatogram(outline = TRUE, fillOutline='grey95', organism = "human", sex = anatogram, fill = 'value') +
      ggplot2::coord_fixed() +
      scale_fill_ddh_c(palette = "gene", breaks = break_points, name = NULL) +
      ggplot2::guides(fill = ggplot2::guide_colorsteps(show.limits = TRUE)) +
      ggplot2::theme_void(base_size = 16, base_family = "Chivo") +
      NULL

    if(card == TRUE){
      body_plot <-
        body_plot +
        ggplot2::labs(x = "") + #, title = "Tissue Distribution", caption = "more ...") +
        NULL
    }

    return(body_plot)
  }
  #error handling
  tryCatch(make_female_anatogram_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

#' Male Anatogram Plot
#'
#' This is a plot function that takes a gene name and returns a male anatogram plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#'
#' @return If no error, then returns a male anatogram plot. If an error is thrown, then will return a bomb plot.
#'
#' @import gganatogram
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_male_anatogram(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_male_anatogram(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_male_anatogram(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_male_anatogram <- function(anatogram = "male",
                                input = list(),
                                card = FALSE){

  male_anatogram <- make_female_anatogram(anatogram = anatogram,
                                          input = input,
                                          card = card)
  return(male_anatogram)
}

##TISSUE GEOM_COL ------------------------------------------
#' Tissue Plot
#'
#' \code{make_tissue} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a tissue plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a tissue plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_tissue(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_tissue(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_tissue(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_tissue <- function(data_gene_tissue = gene_tissue,
                        input = list(),
                        card = FALSE) {
  make_tissue_raw <- function() {
    plot_data <-
      data_gene_tissue %>%
      dplyr::filter_all(dplyr::any_vars(gene_name %in% input$content)) %>%
      dplyr::group_by(organ) %>%
      dplyr::mutate(sum_value = sum(value),
                    organ = stringr::str_replace_all(organ, "_", " "),
                    organ = stringr::str_to_title(organ)) %>%
      dplyr::arrange(dplyr::desc(sum_value))

    if(nrow(plot_data) == 0){return(NULL)}

    if(card == TRUE){
      #get distinct
      top_organs <- plot_data %>%
        dplyr::distinct(organ) %>%
        dplyr::ungroup() %>%
        dplyr::slice_head(n = 10) %>%
        dplyr::pull()
      plot_data <-
        plot_data %>%
        dplyr::filter(organ %in% top_organs)
    }

    plot_draft <-
      ggplot2::ggplot(plot_data,
                      ggplot2::aes(x = value,
                                   y = forcats::fct_reorder(organ, sum_value))) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::scale_x_continuous(expand = c(0, 0), sec.axis = ggplot2::dup_axis()) +
      ggplot2::scale_y_discrete(expand = c(.01, .01)) +
      theme_ddh() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 0, family = "Roboto Slab"),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.line.y = ggplot2::element_blank()) +
      ggplot2::labs(y = NULL)

    ## Version with color mapped to count for single gene queries
    if(length(input$content) == 1){
      plot_complete <-
        plot_draft +
        ggplot2::geom_col(ggplot2::aes(fill = value), width = .82) +
        scale_fill_ddh_c(palette = "gene", guide = "none") +
        ggplot2::labs(x = paste0(stringr::str_c(input$content, collapse = ", "), " Normalized Expression")) +
        ggplot2::theme(plot.margin = ggplot2::margin(0, 15, 0, 0)) # add some space to avoid cutting of labels
    } else {
      plot_complete <-
        plot_draft +
        ggplot2::geom_col(ggplot2::aes(fill = gene_name), width = .82) +
        scale_fill_ddh_d(palette = "gene", shuffle = TRUE, seed = 5L) +
        ggplot2::labs(x = "Sum of Normalized Expression",
                      fill = "Query\nGene") +
        ggplot2::theme(legend.justification = "top")
    }

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        ggplot2::labs(x = "") +
        ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                       legend.position='none') +
        NULL
    }

    return(plot_complete)
  }
  #error handling
  tryCatch(make_tissue_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## EXPRESSION PLOT --------------------------------------------------------
#' Cell Expression Plot
#'
#' \code{make_cellexpression} returns an image of ...
#'
#' Each point shows the ranked expression value across CCLE cell lines.
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a cellexpression plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_cellexpression(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_cellexpression(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_cellexpression(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_cellexpression <- function(data_universal_expression_long = universal_expression_long,
                                data_cell_expression_names = cell_expression_names,
                                data_universal_stats_summary = universal_stats_summary,
                                input = list(),
                                var = "gene",
                                card = FALSE) {
  make_cellexpression_raw <- function() {
    if (var == "gene") {
      plot_initial <-
        data_universal_expression_long %>%
        dplyr::select(dplyr::any_of(c("X1", "gene", "gene_expression"))) %>%
        dplyr::rename("expression_var" = "gene_expression")
      mean <- get_stats(data_set = "expression_gene", var = "mean")
      upper_limit <- get_stats(data_set = "expression_gene", var = "upper")
      lower_limit <- get_stats(data_set = "expression_gene", var = "lower")
      color_type <- "gene"
    } else if (var == "protein") {
      plot_initial <-
        data_universal_expression_long %>%
        dplyr::select(dplyr::any_of(c("X1", "gene", "protein_expression"))) %>%
        dplyr::rename("expression_var" = "protein_expression")
      mean <- get_stats(data_set = "expression_protein", var = "mean")
      upper_limit <- get_stats(data_set = "expression_protein", var = "upper")
      lower_limit <- get_stats(data_set = "expression_protein", var = "lower")
      color_type <- "protein"
    } else {
      stop("declare your variable")
    }

    if (input$type == "gene") {
      plot_data <-
        plot_initial %>%
        dplyr::filter(gene %in% input$content,
                      !is.na(expression_var)) %>%
        dplyr::left_join(data_cell_expression_names, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::select(cell_line, lineage, lineage_subtype, dplyr::everything()) %>%
        dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
        dplyr::mutate(gene_fct = forcats::fct_inorder(gene)) %>%
        ggplot2::ggplot(ggplot2::aes(y = gene_fct,
                                     x = expression_var,
                                     text = paste0("Cell Line: ", cell_line),
                                     color = gene
        ))
    } else if (input$type == "cell") {
      color_type <- "cell"
      plot_data <-
        plot_initial %>%
        dplyr::left_join(data_cell_expression_names, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::select(cell_line, lineage, lineage_subtype, dplyr::everything()) %>%
        dplyr::filter(cell_line %in% input$content,
                      !is.na(expression_var)) %>%
        dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
        dplyr::mutate(cell_fct = forcats::fct_inorder(cell_line)) %>%
        ggplot2::ggplot(ggplot2::aes(y = cell_fct,
                                     x = expression_var,
                                     text = paste0("Gene: ", gene),
                                     color = cell_line))
      mean <- mean_virtual_achilles_cell_line
      upper_limit <- NULL
      lower_limit <- NULL
    }

    plot_complete <-
      plot_data +
      ggplot2::geom_point(alpha = 0.1, shape = "|", size = 12) +
      ggplot2::geom_vline(xintercept = mean, color = "lightgray", linetype = "dashed") +#3SD
      ggplot2::geom_vline(xintercept = lower_limit, color = "lightgray", linetype = "dashed") +#3SD
      ggplot2::geom_vline(xintercept = upper_limit, color = "lightgray", linetype = "dashed") +#3SD
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.01)) +
      ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = 1 / length(input$content)), na.translate = FALSE) +
      ddh::scale_color_ddh_d(palette = color_type) + #req'd to get protein color
      ddh::theme_ddh() +
      ggplot2::theme(
        text = ggplot2::element_text(family = "Nunito Sans"),
        axis.text = ggplot2::element_text(family = "Roboto Slab"),
        axis.text.y = ggplot2::element_text(size = 18),
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        legend.position = "none"
      ) +
      NULL

    if(length(input$content) == 1){
      plot_complete  <-
        plot_complete +
        ggplot2::theme(axis.text.y = ggplot2::element_blank()) +
        ggplot2::labs(x = paste0(stringr::str_c(input$content, collapse = ", "), " ", var, " levels"))
    } else {
      plot_complete <-
        plot_complete +
        ggplot2::labs(x = paste0(var, " levels"), color = "Query")
    }

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        ggplot2::labs(x = "")
    }
    return(plot_complete)
  }
  #error handling
  tryCatch(print(make_cellexpression_raw()),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

# G-EXPvP-EXP ----------------------------------
#' Gene Expression versus Protein Expression
#'
#' Each point shows the gene expression value compared to the protein expression value for gene within a given cell line. The Pearson correlation coefficient and the p-values are provided in the top-left corner of the plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_cellgeneprotein(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_cellgeneprotein(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_cellgeneprotein(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_cellgeneprotein <- function(data_universal_expression_long = universal_expression_long,
                                 data_cell_expression_names = cell_expression_names,
                                 input = list(),
                                 card = FALSE) {
  make_cellgeneprotein_raw <- function() {
    if (input$type == "gene") {
      plot_initial <-
        data_universal_expression_long %>%
        dplyr::filter(gene %in% input$content,
                      !is.na(gene_expression),
                      !is.na(protein_expression)) %>%
        dplyr::left_join(data_cell_expression_names, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::select(cell_line, lineage, lineage_subtype, dplyr::everything()) %>%
        dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
        ggplot2::ggplot(ggplot2::aes(x = gene_expression,
                                     y = protein_expression,
                                     text = paste0("Cell Line: ", cell_line),
                                     color = gene,
                                     group = gene)
        )
    } else if (input$type == "cell") {
      plot_initial <-
        data_universal_expression_long %>%
        dplyr::left_join(data_cell_expression_names, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::select(cell_line, lineage, lineage_subtype, dplyr::everything()) %>%
        dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
        dplyr::filter(cell_line %in% input$content,
                      !is.na(gene_expression),
                      !is.na(protein_expression)) %>%
        ggplot2::ggplot(ggplot2::aes(x = gene_expression,
                                     y = protein_expression,
                                     text = paste0("Gene: ", gene),
                                     color = cell_line,
                                     group = cell_line)
        )
    }

    plot_complete <-
      plot_initial +
      ggplot2::geom_point(alpha = 0.4) +
      #add geom to drop linear regression line?
      ggplot2::geom_hline(yintercept = 0, color = "lightgray") +
      ggplot2::geom_vline(xintercept = 0, color = "lightgray") +
      # smooth line
      ggplot2::geom_smooth(method = "lm",
                           se = TRUE) +
      # R coefs
      {if(card == FALSE)ggpubr::stat_cor(digits = 3)} +
      scale_color_ddh_d(palette = input$type) +
      theme_ddh() +
      ggplot2::theme(
        text = ggplot2::element_text(family = "Nunito Sans"),
        axis.text = ggplot2::element_text(family = "Roboto Slab")
      ) +
      NULL

    if(length(input$content) == 1){
      plot_complete  <-
        plot_complete +
        ggplot2::labs(x = paste0(input$content, " Gene Expression"), y = paste0(input$content, " Protein Expression")) +
        ggplot2::theme(legend.position = "none")
    } else if(input$type == "gene") {
      plot_complete <-
        plot_complete +
        ggplot2::labs(x = "Gene Expression", y = "Protein Expression", color = "Query \nGene")
    } else if(input$type == "cell") {
      plot_complete <-
        plot_complete +
        ggplot2::labs(x = "Gene Expression", y = "Protein Expression", color = "Query \nCell Line")
    }

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        ggplot2::labs(x = "", y = "") +
        ggplot2::theme(legend.position='none')
    }
    return(plot_complete)
  }
  #error handling
  tryCatch(print(make_cellgeneprotein_raw()),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## CELL DEPS --------------------------------------------------------------------
#' Dependency Curve Plot
#'
#' Each point shows the ranked dependency score ordered from low to high scores. Dependency scores less than -1 indicate a gene that is essential within a cell line. Dependency scores close to 0 mean no changes in fitness when the gene is knocked out. Dependency scores greater than 1 indicate gene knockouts lead to a gain in fitness.
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a celldeps plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_celldeps(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_celldeps(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_celldeps(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_celldeps <- function(data_universal_achilles_long = universal_achilles_long,
                          data_universal_prism_long = universal_prism_long,
                          data_universal_stats_summary = universal_stats_summary,
                          data_cell_expression_names = cell_expression_names,
                          input = list(),
                          card = FALSE,
                          lineplot = FALSE,
                          scale = NULL) {#scale is expecting 0 to 1
  make_celldeps_raw <- function() {
    if(input$type == "gene") {
      aes_var <- rlang::sym("name")
      var_title <- "Gene"
      ylab <- "Dependecy Score"
      mean <- get_stats(data_set = "achilles", var = "mean")

      plot_data <-
        data_universal_achilles_long %>% #plot setup
        dplyr::filter(gene %in% input$content) %>%
        dplyr::left_join(data_cell_expression_names, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::group_by(gene) %>%
        dplyr::arrange(dep_score) %>%
        dplyr::mutate(
          rank = 1:dplyr::n(),
          med = median(dep_score, na.rm= TRUE)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::rename(name = gene)

    } else if(input$type == "compound") {
      aes_var <- rlang::sym("name")
      var_title <- "Compound"
      ylab <- "Log2FC"
      mean <- get_stats(data_set = "prism", var = "mean")

      plot_data <-
        data_universal_prism_long %>% #plot setup
        dplyr::filter(name %in% input$content) %>%
        dplyr::left_join(data_cell_expression_names, by = c("x1" = "X1")) %>%
        dplyr::select(-1) %>%
        dplyr::group_by(name) %>%
        dplyr::arrange(log2fc) %>%
        dplyr:: mutate(
          rank = 1:n(),
          med = median(log2fc, na.rm= TRUE)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::rename(dep_score = log2fc) #rename for graph

    } else { #cell lines
      aes_var <- rlang::sym("cell_line")
      var_title <- "Gene"
      ylab <- "Dependecy Score"
      mean <- get_stats(data_set = "achilles_cell", var = "mean")

      plot_data <-
        data_universal_achilles_long %>% #plot setup
        dplyr::left_join(data_cell_expression_names, by = "X1") %>%
        dplyr::filter(cell_line %in% input$content) %>%
        dplyr::select(-X1) %>%
        dplyr::group_by(cell_line) %>%
        dplyr::arrange(dep_score) %>%
        dplyr:: mutate(
          rank = 1:dplyr::n(),
          med = median(dep_score, na.rm= TRUE)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::rename(name = gene)
    }

    if(!is.null(scale) & !card){
      plot_data <-
        plot_data %>%
        dplyr::slice_sample(prop = scale)
    }

    if(card) {
      if(is.null(scale)) {
        scale <- 0.3
      }
      if(input$type == "gene") {
        plot_data <-
          plot_data %>%
          dplyr::group_by(name) %>%
          dplyr::sample_n(scale*dplyr::n()) %>%
          dplyr::ungroup()
      } else {
        plot_data <-
          plot_data %>%
          dplyr::group_by(cell_line) %>%
          dplyr::sample_n(scale*dplyr::n()) %>%
          dplyr::ungroup()
      }
    }

    plot_complete <-
      plot_data %>%
      ggplot2::ggplot(ggplot2::aes(x = rank,
                                   y = dep_score,
                                   text = glue::glue('{var_title}: {name}\nCell Line: {cell_line}'),
                                   color = forcats::fct_reorder(!!aes_var, med),
                                   fill = forcats::fct_reorder(!!aes_var, med)
      )) +
      ## dot/line plot
      {if(!card & !lineplot)ggplot2::geom_point(size = 1.1, stroke = .25, alpha = 0.6)} +
      {if(input$type == "gene" & (card | lineplot))ggplot2::geom_line(ggplot2::aes(group = name))} +
      {if(input$type == "cell" & (card | lineplot))ggplot2::geom_line(ggplot2::aes(group = cell_line))} +
      ## indicator lines dep. score
      ggplot2::geom_hline(yintercept = mean, linetype = "dashed", color = "grey80") +
      ggplot2::geom_hline(yintercept = 1, size = .2, color = "grey70") +
      ggplot2::geom_hline(yintercept = -1, size = .2, color = "grey70") +
      ## scales + legends
      #scale_x_discrete(expand = expansion(mult = 0.02), na.translate = FALSE) +
      scale_color_ddh_d(palette = input$type) +
      scale_fill_ddh_d(palette = input$type) +
      ggplot2::guides(
        color = ggplot2::guide_legend(reverse = TRUE, override.aes = list(size = 4, stroke = .8)),
        fill = ggplot2::guide_legend(reverse = TRUE, override.aes = list(size = 3.8, stroke = .8))
      ) +
      ## titles
      ggplot2::labs(
        x = NULL,
        y = ylab,
        color = "Query",
        fill = "Query"
      ) +
      ## theme changes
      theme_ddh() +
      ggplot2::theme(
        text = ggplot2::element_text(family = "Nunito Sans"),
        axis.text = ggplot2::element_text(family = "Roboto Slab"),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank()
      ) +
      NULL

    ##only plot legend in case of more than 1 gene selected
    if(length(input$content) == 1){
      plot_complete  <-
        plot_complete +
        ggplot2::theme(legend.position = "none")
      plot_complete
    } else {
      plot_complete
    }

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        ggplot2::labs(x = "") +
        ggplot2::theme(legend.position = "none")
    }

    return(plot_complete)
  }
  #error handling
  tryCatch(make_celldeps_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## BAR PLOT --------------------------------------------------------------------
#' Cell Dependencies Bar Plot
#'
#' Each bar shows the dependency scores of the queried genes in a cell line. Dependency scores less than -1 indicate a gene that is essential within a cell line. Dependency scores close to 0 mean no changes in fitness when the gene is knocked out. Dependency scores greater than 1 indicate gene knockouts lead to a gain in fitness.
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a cellbar plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_cellbar(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_cellbar(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_cellbar(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_cellbar <- function(data_universal_achilles_long = universal_achilles_long,
                         data_universal_prism_long = universal_prism_long,
                         data_cell_expression_names = cell_expression_names,
                         data_universal_stats_summary = universal_stats_summary,
                         input = list(),
                         card = FALSE,
                         scale = NULL) {
  make_cellbar_raw <- function() {
    if(input$type == "gene") {
      aes_var <- rlang::sym("name")
      var_title <- "Gene"
      ylab <- "Dependecy Score"
      mean <- get_stats(data_set = "achilles", var = "mean")

      plot_data <-
        data_universal_achilles_long %>% #plot setup
        dplyr::filter(gene %in% input$content) %>%
        dplyr::left_join(data_cell_expression_names, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::group_by(gene) %>%
        dplyr::arrange(dep_score) %>%
        dplyr::mutate(
          rank = 1:dplyr::n(),
          med = median(dep_score, na.rm= TRUE)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::rename(name = gene) %>%
        dplyr::mutate(rank = as.integer(forcats::fct_reorder(cell_line, dep_score)))

    } else if(input$type == "compound") {
      aes_var <- rlang::sym("name")
      var_title <- "Compound"
      ylab <- "Log2FC"
      mean <- get_stats(data_set = "prism", var = "mean")

      plot_data <-
        data_universal_prism_long %>% #plot setup
        dplyr::filter(name %in% input$content) %>%
        dplyr::left_join(data_cell_expression_names, by = c("x1" = "X1")) %>%
        dplyr::select(-1) %>%
        dplyr::group_by(name) %>%
        dplyr::arrange(log2fc) %>%
        dplyr::mutate(
          rank = 1:dplyr::n(),
          med = median(log2fc, na.rm= TRUE)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::rename(dep_score = log2fc) %>% #rename for graph
        dplyr::mutate(rank = as.integer(forcats::fct_reorder(name, dep_score)))

    } else { #cell lines
      aes_var <- rlang::sym("cell_line")
      var_title <- "Gene"
      ylab <- "Dependecy Score"
      mean <- get_stats(data_set = "achilles_cell", var = "mean")

      plot_data <-
        data_universal_achilles_long %>% #plot setup
        dplyr::left_join(data_cell_expression_names, by = "X1") %>%
        dplyr::filter(cell_line %in% input$content) %>%
        dplyr::select(-X1) %>%
        dplyr::group_by(cell_line) %>%
        dplyr::arrange(dep_score) %>%
        dplyr::mutate(
          rank = 1:dplyr::n(),
          med = median(dep_score, na.rm= TRUE)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::rename(name = gene) %>%
        dplyr::mutate(rank = as.integer(forcats::fct_reorder(name, dep_score)))

    }

    if(card) {
      if(is.null(scale)) {
        scale <- 0.3
      }
      if(input$type == "gene") {
        plot_data <-
          plot_data %>%
          dplyr::group_by(name) %>%
          dplyr::sample_n(scale*dplyr::n()) %>%
          dplyr::ungroup()
      } else {
        plot_data <-
          plot_data %>%
          dplyr::group_by(cell_line) %>%
          dplyr::sample_n(scale*dplyr::n()) %>%
          dplyr::ungroup()
      }
    }

    plot_complete <-
      plot_data %>%
      ggplot2::ggplot(ggplot2::aes(x = rank,
                                   y = dep_score,
                                   text = glue::glue('{var_title}: {name}\nCell Line: {cell_line}'),
                                   color = forcats::fct_reorder(!!aes_var, med),
                                   fill = forcats::fct_reorder(!!aes_var, med)
      )) +
      ## bar plot
      ggplot2::geom_bar(stat = "identity", width = 0.5) +
      ## indicator lines dep. score
      ggplot2::geom_hline(yintercept = mean, linetype = "dashed", color = "grey80") +
      ggplot2::geom_hline(yintercept = 1, size = .2, color = "grey70") +
      ggplot2::geom_hline(yintercept = -1, size = .2, color = "grey70") +
      ggplot2::geom_hline(yintercept = 0, size = .2, color = "grey70") +
      ## scales + legends
      ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = 0.02), na.translate = FALSE) +
      scale_color_ddh_d(palette = input$type) +
      scale_fill_ddh_d(palette = input$type) +
      ggplot2::guides(
        color = ggplot2::guide_legend(reverse = TRUE, override.aes = list(size = 4, stroke = .8)),
        fill = ggplot2::guide_legend(reverse = TRUE, override.aes = list(size = 3.8, stroke = .8))
      ) +
      ## titles
      ggplot2::labs(
        x = NULL,
        y = ylab,
        color = "Query",
        fill = "Query"
      ) +
      ## theme changes
      theme_ddh() + #base_size = 15 default
      ggplot2::theme(
        text = ggplot2::element_text(family = "Nunito Sans"),
        axis.text = ggplot2::element_text(family = "Roboto Slab"),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank()
      ) +
      NULL

    ## only plot legend in case of more than 1 gene selected
    if(length(input$content) == 1){
      plot_complete  <-
        plot_complete +
        ggplot2::theme(legend.position = "none")
    }

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        ggplot2::labs(x = "") +
        ggplot2::theme(legend.position = "none")
    }
    return(plot_complete)
  }
  #error handling
  tryCatch(make_cellbar_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## DENSITY PLOT ----------------------------------------------------------------
#' Cell Dependencies Density Plot
#'
#' Kernel density estimate of dependency scores. Dependency scores across all cell lines for queried genes, revealing overall influence of a gene on cellular fitness. The interval indicates the 95% quantile of the data, the dot indicates the median dependency score. The gray background highlights weak dependency values between -1 and 1.
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a cellbins plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#'
#' make_cellbins(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_cellbins(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_cellbins(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_cellbins <- function(data_universal_achilles_long = universal_achilles_long,
                          data_universal_prism_long = universal_prism_long,
                          data_cell_expression_names = cell_expression_names,
                          input = list(),
                          card = FALSE) {
  make_cellbins_raw <- function() {
    if(input$type == "gene") {
      xlab <- "Dependency Score (distribution)"
      aes_var <- rlang::sym("name")

      plot_data <-
        data_universal_achilles_long %>% #plot setup
        dplyr::filter(gene %in% input$content) %>%
        dplyr::left_join(data_cell_expression_names, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::group_by(gene) %>%
        dplyr::arrange(dep_score) %>%
        dplyr::mutate(med = median(dep_score, na.rm= TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(!is.na(dep_score)) %>%
        dplyr::rename(name = gene)

    } else if(input$type == "compound") {
      xlab <- "Log2FC (distribution)"
      aes_var <- rlang::sym("name")

      plot_data <-
        data_universal_prism_long %>% #plot setup
        dplyr::filter(name %in% input$content) %>%
        dplyr::left_join(data_cell_expression_names, by = c("x1" = "X1")) %>%
        dplyr::select(-1) %>%
        dplyr::group_by(name) %>%
        dplyr::arrange(log2fc) %>%
        dplyr::mutate(
          rank = 1:dplyr::n(),
          med = median(log2fc, na.rm= TRUE)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::rename(dep_score = log2fc) #rename for graph

    } else if(input$type == "cell") {
      xlab <- "Dependency Scores (distribution)"
      aes_var <- rlang::sym("cell_line")

      plot_data <-
        data_universal_achilles_long %>% #plot setup
        dplyr::left_join(data_cell_expression_names, by = "X1") %>%
        dplyr::filter(cell_line %in% input$content) %>%
        dplyr::select(-X1) %>%
        dplyr::group_by(cell_line) %>%
        dplyr::arrange(dep_score) %>%
        dplyr::mutate(med = median(dep_score, na.rm= TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(!is.na(dep_score)) %>%
        dplyr::rename(name = gene)

    } else {
      return("stop! delcare your type")
    }

    colors <- ddh_pal_d(palette = input$type)(length(input$content))

    plot_complete <-
      plot_data %>%
      dplyr::filter(!is.na(dep_score)) %>%
      ggplot2::ggplot() +
      ## annotation range -1 to 1
      ggplot2::geom_rect(
        xmin = -1, xmax = 1,
        ymin = -Inf, ymax = Inf,
        fill = "grey95",
        show.legend = FALSE
      ) +
      ## indicator line y axis
      ggplot2::geom_linerange(
        ggplot2::aes(xmin = -Inf, xmax = med,
                     y = forcats::fct_reorder(!!aes_var, -med),
                     color = med < -1),
        linetype = "dotted",
        size = .2,
        show.legend = FALSE
      ) +
      ## density curves via {ggdist}
      ggdist::stat_halfeye(ggplot2::aes(x = dep_score,
                                        y = forcats::fct_reorder(!!aes_var, -med),
                                        fill = stat(abs(x) > 1),
                                        point_fill = ggplot2::after_scale(fill)),
                           .width = c(.025, .975),
                           color = "black",
                           shape = 21,
                           stroke = .7,
                           point_size = 2) +
      ## zero line
      ggplot2::geom_vline(
        xintercept = 0,
        color = "grey80",
        linetype = "dashed",
        show.legend = FALSE
      ) +
      ## titles
      ggplot2::labs(
        x = xlab,
        y = NULL,
        color = "Query",
        fill = "Query"
      ) +
      ## scales + legends
      ggplot2::scale_y_discrete(expand = c(.03, .03)) +
      ggplot2::scale_color_manual(values = c("grey70", colors)) +
      ggplot2::scale_fill_manual(values = c("grey70", colors)) +
      ggplot2::guides(
        color = ggplot2::guide_legend(size = 1, reverse = TRUE),
        fill = ggplot2::guide_legend(size = 1, reverse = TRUE)
      ) +
      ## theme changes
      theme_ddh() +
      ggplot2::theme(
        text = ggplot2::element_text(family = "Nunito Sans"),
        legend.position = "none",
        axis.line.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(family = "Roboto Slab", size = 18),
        axis.text.x = ggplot2::element_text(size = 12, color = "grey30"),
        axis.title = ggplot2::element_text(size = 15)
      ) +
      NULL

    if(card) {
      plot_complete <-
        plot_complete +
        ggplot2::labs(x = "") + #too long of a xlabel
        NULL
    }
    return(plot_complete)
  }
  #error handling
  tryCatch(print(make_cellbins_raw()), #add print to catch ggplot.print errors
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## LINEAGE LINERANGE PLOT ------------------------------------------------------
#' Dependency Lineage Plot
#'
#' Each point shows the mean dependency score for the gene query within a given cell lineage. The intervals show the 5% quantiles centered on the median, the interquartile ranges, and the 95% quantiles. The gray background highlights weak dependency values between -1 and 1.
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a lineage plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_lineage(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_lineage(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), highlight = TRUE)
#' make_lineage(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_lineage(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_lineage <- function(data_universal_achilles_long = universal_achilles_long,
                         data_universal_prism_long = universal_prism_long,
                         data_cell_expression_names = cell_expression_names,
                         input = list(),
                         card = FALSE,
                         highlight = FALSE) {
  make_lineage_raw <- function() {
    if(input$type == "gene" ) {
      xlab <- "Dependency Score"
      title_var <- glue::glue('Cell lineage dependencies for {stringr::str_c(input$content, collapse = ", ")}')

      data_full <-
        data_universal_achilles_long %>% #plot setup
        dplyr::filter(gene %in% input$content) %>%
        dplyr::left_join(data_cell_expression_names, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::mutate_at("lineage", function(str) {
          str <- stringr::str_replace_all(str, "\\_", " ")
          str <- stringr::str_to_title(str)
          return(str)
        }) %>%
        tidyr::drop_na(lineage) %>%
        tidyr::drop_na(dep_score) %>%
        dplyr::group_by(lineage) %>%
        dplyr::mutate(mean = mean(dep_score)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(lineage = forcats::fct_reorder(lineage, -mean))

      data_mean <- data_full %>%
        dplyr::group_by(lineage) %>%
        dplyr::summarize(dep_score = mean(dep_score))

    } else if(input$type == "compound") {
      xlab <- "Log2FC"
      title_var <- glue::glue('Cell lineage dependencies for {stringr::str_c(input$content, collapse = ", ")}')

      data_full <-
        data_universal_prism_long %>% #plot setup
        dplyr::filter(name %in% input$content) %>%
        dplyr::left_join(data_cell_expression_names, by = c("x1" = "X1")) %>%
        dplyr::select(-1) %>%
        dplyr::mutate_at("lineage", function(str) {
          str <- stringr::str_replace_all(str, "\\_", " ")
          str <- stringr::str_to_title(str)
          return(str)
        }) %>%
        tidyr::drop_na(lineage) %>%
        tidyr::drop_na(log2fc) %>%
        dplyr::group_by(lineage) %>%
        dplyr::mutate(mean = mean(log2fc)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(lineage = forcats::fct_reorder(lineage, -mean)) %>%
        dplyr::rename(dep_score = log2fc) #rename for graph

      data_mean <-
        data_full %>%
        dplyr::group_by(lineage) %>%
        dplyr::summarize(dep_score = mean(dep_score))

    } else {
      stop("declare your type") }

    if(nrow(data_full) == 0){return(NULL)}

    if(highlight) {
      stats_data <- data_full %>%
        dplyr::group_by(lineage) %>%
        dplyr::filter(dplyr::n() > 1) %>%
        dplyr::ungroup() %>%
        ggpubr::compare_means(dep_score ~ lineage, data = .,
                              ref.group = ".all.", method = "t.test") %>%
        dplyr::filter(p.adj < 0.05)
    }

    if(card) {
      most_negative <-
        data_mean %>%
        dplyr::slice_min(dep_score, n = 6) %>%
        dplyr::pull(lineage)
    }

    plot_complete <-
      data_mean %>%
      {if(card)dplyr::filter(., lineage %in% most_negative) else .} %>%
      ggplot2::ggplot(ggplot2::aes(dep_score, lineage)) +
      ## annotation range -1 to 1
      # geom_rect(
      #   xmin = -1, xmax = 1,
      #   ymin = -Inf, ymax = Inf,
      #   fill = "grey95"
      # ) +
      ## zero line
      ggplot2::geom_vline(
        xintercept = 0,
        color = "grey80",
        linetype = "dashed"
      ) +
      ## indicator lines lineages
      ggplot2::geom_linerange(
        ggplot2::aes(xmin = -Inf, xmax = dep_score),
        color = "grey60",
        linetype = "dotted"
      ) +
      ## lineranges as "boxplots"
      ggdist::stat_interval(
        data = data_full %>% {if(card)dplyr::filter(., lineage %in% most_negative) else .},
        orientation = "horizontal",
        .width = c(.05, .5, .95)
      ) +
      ## dot indicating mean
      ggplot2::geom_point(
        color = "black", fill = "white",
        shape = 21, stroke = .5,
        size = 1.8
      ) +
      ## scales + legends
      ggplot2::scale_x_continuous(
        sec.axis = ggplot2::dup_axis()
      ) +
      scale_color_ddh_d(
        palette = input$type,
        shuffle = TRUE, seed = 5L, ## to return "correctly" ordered, sequential colors
        labels = c("95%", "50%", "5%"),#of the data fall in these ranges
        name = ""
      ) +
      ggplot2::guides(color = ggplot2::guide_legend(reverse = TRUE)) +
      ## titles
      ggplot2::labs(
        x = xlab,
        y = NULL #,
        #title = title_var
      ) +
      {if(highlight)gghighlight::gghighlight(lineage %in% stats_data$group2, use_direct_label = FALSE)} + # toggle
      ## theme changes
      theme_ddh(grid = "none") +
      ggplot2::theme(
        legend.position = "top",
        axis.line.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(family = "Roboto Slab"),
        axis.text.x = ggplot2::element_text(size = 12, color = "grey30"),
        axis.title = ggplot2::element_text(size = 15),
        axis.title.x.bottom = ggplot2::element_blank(),
        plot.title.position = "plot"
      ) +
      NULL

    if(card) {
      plot_complete <-
        plot_complete +
        ggplot2::scale_y_discrete(labels = scales::label_wrap(10)) + #to prevent long lines and squished plots
        ggplot2::scale_x_continuous(breaks = scales::breaks_extended(n = 3)) +
        ggplot2::theme(plot.title = ggplot2::element_blank())
    }
    return(plot_complete)
  }
  #error handling
  tryCatch(make_lineage_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## SUBLINE RANGE PLOT ---------------------------------------------------
#' Cell Dependency Sublineage Plot
#'
#' Each point shows the mean dependency score for the gene query within a given cell lineage. The intervals show the 5% quantiles centered on the median, the interquartile ranges, and the 95% quantiles. The gray background highlights weak dependency values between -1 and 1.
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a sublineage plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_sublineage(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_sublineage(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), highlight = TRUE)
#' make_sublineage(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_sublineage(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_sublineage <- function(data_universal_achilles_long = universal_achilles_long,
                            data_universal_prism_long = universal_prism_long,
                            data_cell_expression_names = cell_expression_names,
                            input = list(),
                            card = FALSE,
                            highlight = FALSE) {
  make_sublineage_raw <- function() {
    if(input$type == "gene") {
      xlab <- "Dependency Score"
      title_var <- glue::glue('Cell sub-lineage dependencies for {stringr::str_c(input$content, collapse = ", ")}')

      data_full <-
        data_universal_achilles_long %>% #plot setup
        dplyr::filter(gene %in% input$content) %>%
        dplyr::left_join(data_cell_expression_names, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::mutate_at("lineage_subtype", function(str) {
          str <- stringr::str_replace_all(str, "\\_", " ")
          str <- dplyr::if_else(stringr::str_detect(str, "^[:lower:]"), stringr::str_to_title(str), str)
          return(str)
        })  %>%
        tidyr::drop_na(lineage_subtype) %>%
        tidyr::drop_na(dep_score) %>%
        dplyr::group_by(lineage_subtype) %>%
        dplyr::mutate(mean = mean(dep_score)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(lineage_subtype = forcats::fct_reorder(lineage_subtype, -mean))

      data_mean <- data_full %>%
        dplyr::group_by(lineage_subtype) %>%
        dplyr::summarize(dep_score = mean(dep_score))

    } else if(input$type == "compound") {
      xlab <- "Log2FC"
      title_var <- glue::glue('Cell sub-lineage dependencies for {stringr::str_c(input$content, collapse = ", ")}')

      data_full <-
        data_universal_prism_long %>% #plot setup
        dplyr::filter(name %in% input$content) %>%
        dplyr::left_join(data_cell_expression_names, by = c("x1" = "X1")) %>%
        dplyr::select(-1) %>%
        dplyr::mutate_at("lineage_subtype", function(str) {
          str <- stringr::str_replace_all(str, "\\_", " ")
          str <- dplyr::if_else(stringr::str_detect(str, "^[:lower:]"), stringr::str_to_title(str), str)
          return(str)
        })  %>%
        tidyr::drop_na(lineage_subtype) %>%
        tidyr::drop_na(log2fc) %>%
        dplyr::group_by(lineage_subtype) %>%
        dplyr::mutate(mean = mean(log2fc)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(lineage_subtype = forcats::fct_reorder(lineage_subtype, -mean)) %>%
        dplyr::rename(dep_score = log2fc) #rename for graph

      data_mean <- data_full %>%
        dplyr::group_by(lineage_subtype) %>%
        dplyr::summarize(dep_score = mean(dep_score))

    } else {
      stop("declare your type") }

    if(nrow(data_full) == 0){return(NULL)}

    if(highlight) {
      stats_data <- data_full %>%
        dplyr::group_by(lineage_subtype) %>%
        dplyr::filter(dplyr::n() > 1) %>%
        dplyr::ungroup() %>%
        ggpubr::compare_means(dep_score ~ lineage_subtype, data = .,
                              ref.group = ".all.", method = "t.test") %>%
        dplyr::filter(p.adj < 0.05)
    }

    if(card) {
      most_negative <-
        data_mean %>%
        dplyr::slice_min(dep_score, n = 6) %>%
        dplyr::pull(lineage_subtype)
    }

    plot_complete <-
      data_mean %>%
      {if(card)dplyr::filter(., lineage_subtype %in% most_negative) else .} %>%
      ggplot2::ggplot(ggplot2::aes(dep_score, lineage_subtype)) +
      ## annotation range -1 to 1
      # geom_rect(
      #   xmin = -1, xmax = 1,
      #   ymin = -Inf, ymax = Inf,
      #   fill = "grey95"
      # ) +
      ## zero line
      ggplot2::geom_vline(
        xintercept = 0,
        color = "grey80",
        linetype = "dashed"
      ) +
      ## indicator lines sublineages
      ggplot2::geom_linerange(
        ggplot2::aes(xmin = -Inf, xmax = dep_score),
        color = "grey60",
        linetype = "dotted"
      ) +
      ## lineranges as "boxplots"
      ggdist::stat_interval(
        data = data_full %>% {if(card)dplyr::filter(., lineage_subtype %in% most_negative) else .},
        orientation = "horizontal",
        .width = c(.05, .5, .95)
      ) +
      ## dot indicating mean
      ggplot2::geom_point(
        color = "black", fill = "white",
        shape = 21, stroke = .5,
        size = 1.8
      ) +
      ## scales + legends
      ggplot2::scale_x_continuous(
        sec.axis = ggplot2::dup_axis()
      ) +
      scale_color_ddh_d(
        palette = input$type,
        shuffle = TRUE, seed = 5L, ## to return "correctly" ordered, sequential colors
        labels = c("95%", "50%", "5%"), #of the data fall in these ranges
        name = ""
      ) +
      ggplot2::guides(color = ggplot2::guide_legend(reverse = TRUE)) +
      ## titles
      ggplot2::labs(
        x = xlab,
        y = NULL #,
        #title = title_var
      ) +
      {if(highlight)gghighlight::gghighlight(lineage_subtype %in% stats_data$group2, use_direct_label = FALSE)} + # toggle
      ## theme changes
      theme_ddh(grid = "none") +
      ggplot2::theme(
        legend.position = "top",
        axis.line.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(family = "Roboto Slab"),
        axis.text.x = ggplot2::element_text(size = 12, color = "grey30"),
        axis.title = ggplot2::element_text(size = 15),
        axis.title.x.bottom = ggplot2::element_blank(),
        plot.title.position = "plot"
      ) +
      NULL

    if(card == TRUE) {
      plot_complete <-
        plot_complete +
        ggplot2::scale_y_discrete(labels = scales::label_wrap(10)) + #to prevent long lines and squished plots
        ggplot2::scale_x_continuous(breaks = scales::breaks_extended(n = 3)) +
        ggplot2::theme(plot.title = ggplot2::element_blank())
    }

    return(plot_complete)
  }
  #error handling
  tryCatch(make_sublineage_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## CORRELATION PLOT FOR CELL DEPS--------------------------------------------------------
#' Co-essentiality Correlation Plot
#'
#' Each point shows the ranked correlation value ordered from high to low for each query. Correlation values outside the solid gray lines indicate the gene has a correlation value greater than the mean.
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a correlation plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_correlation(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_correlation(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_correlation(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_correlation <- function(data_gene_achilles_cor_nest = gene_achilles_cor_nest,
                             #data_achilles_cell_line_cor_nest = achilles_cell_line_cor_nest,
                             data_prism_cor_nest = prism_cor_nest,
                             data_universal_stats_summary = universal_stats_summary,
                             input = list(),
                             card = FALSE,
                             scale = NULL) { #no card option, but need this to prevent error
  make_correlation_raw <- function() {
    if(input$type == "gene") {
      mean <- get_stats(data_set = "achilles", var = "mean")
      upper_limit <- get_stats(data_set = "achilles", var = "upper")
      lower_limit <- get_stats(data_set = "achilles", var = "lower")
      var <- rlang::sym("fav_gene") #from https://rlang.r-lib.org/reference/quasiquotation.html
      label_var <- "Gene Rank"
      text_var <- "Gene"
      content_var <- glue::glue_collapse(input$content, sep = ", ")

      plot_data <-
        data_gene_achilles_cor_nest %>%
        dplyr::filter(fav_gene %in% input$content) %>%
        tidyr::unnest(data) %>%
        dplyr::group_by(fav_gene) %>%
        dplyr::arrange(dplyr::desc(r2)) %>%
        dplyr::mutate(
          rank = 1:dplyr::n(),
          med = median(r2, na.rm= TRUE)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::rename(name = gene) #for plot name

    } else if(input$type == "compound") {
      mean <- get_stats(data_set = "prism", var = "mean")
      upper_limit <- get_stats(data_set = "prism", var = "upper")
      lower_limit <- get_stats(data_set = "prism", var = "lower")
      var <- rlang::sym("fav_drug") #from https://rlang.r-lib.org/reference/quasiquotation.html
      label_var <- "Drug Rank"
      text_var <- "Compound"
      content_var <- glue::glue_collapse(input$content, sep = ", ")

      plot_data <-
        data_prism_cor_nest %>%
        dplyr::filter(fav_drug %in% input$content) %>%
        tidyr::unnest(data) %>%
        dplyr::group_by(fav_drug) %>%
        dplyr::arrange(dplyr::desc(r2)) %>%
        dplyr::mutate(
          rank = 1:dplyr::n(),
          med = median(r2, na.rm= TRUE)
        ) %>%
        dplyr::ungroup()
    } #no else for cell lines

    if(!is.null(scale) & !card){
      plot_data <-
        plot_data %>%
        dplyr::slice_sample(prop = scale)
    }

    if(card) {
      if(is.null(scale)) {
        scale <- 0.3
      }
      plot_data <- plot_data %>%
        dplyr::group_by(name) %>%
        dplyr::sample_n(scale*dplyr::n()) %>% # scale is also used here
        dplyr::ungroup()
    }

    plot_complete <-
      plot_data %>%
      ggplot2::ggplot() +
      ## sd square
      ggplot2::geom_hline(yintercept = upper_limit, color = "gray80", linetype = "dashed") +
      ggplot2::geom_hline(yintercept = mean, color = "gray80", linetype = "dashed") +
      ggplot2::geom_hline(yintercept = lower_limit, color = "gray80", linetype = "dashed") +
      ggplot2::annotate("text", x = Inf, y = 0.005, label = label_var, color = "gray40", hjust = 1.15 ,vjust = 0) +
      ## dot plot
      ggplot2::geom_point(ggplot2::aes(x = rank,
                                       y = r2,
                                       text = glue::glue('{text_var}: {name}'),
                                       color = forcats::fct_reorder(!!var, med), #from https://rlang.r-lib.org/reference/quasiquotation.html
                                       fill = forcats::fct_reorder(!!var, med) #from https://rlang.r-lib.org/reference/quasiquotation.html
      ),
      size = 1.1, stroke = .1, alpha = 0.4) +
      ## scales + legends
      ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = 0.02), na.translate = FALSE) +
      scale_color_ddh_d(palette = input$type) +
      scale_fill_ddh_d(palette = input$type) +
      ggplot2::guides(
        color = ggplot2::guide_legend(reverse = TRUE, override.aes = list(size = 5)),
        fill = ggplot2::guide_legend(reverse = TRUE, override.aes = list(size = 5))
      ) +
      ## labels
      ggplot2::labs(
        x = NULL,
        y = glue::glue("{text_var} correlations with {content_var}"),
        color = glue::glue("Query {text_var}"),
        fill = glue::glue("Query {text_var}")
      ) +
      ## theme changes
      theme_ddh() +
      ggplot2::theme(
        text = ggplot2::element_text(family = "Nunito Sans"),
        axis.text = ggplot2::element_text(family = "Roboto Slab"),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank()
      ) +
      NULL

    ##change axis in case of more than 1 gene selected
    if(length(input$content) == 1){ #fix me
      plot_complete  <-
        plot_complete +
        ggplot2::guides(color = "none",
                        fill = "none")
    }

    return(plot_complete)
  }
  #error handling
  tryCatch(make_correlation_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

## NETWORK PLOT FOR GENE-PATHWAYS --------------------------------------------------------
#' Gene-Pathway co-essentiality plot
#'
#' @param input Expecting a list containing content variable.
#' @return If no error, then returns a network plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_gene_pathways_components_network(input = list(content = 'ROCK1'))
#' \dontrun{
#' make_gene_pathways_components_network(input = list(content = 'ROCK1'))
#' }
make_gene_pathways_components_network <- function(data_universal_achilles_long = universal_achilles_long,
                                                  input = list(),
                                                  highlight = NULL,
                                                  show_labels = FALSE,
                                                  fontsize = 3) {

  make_gene_pathways_components_network_raw <- function() {
    plot_data <- make_gene_pathways_components(input = input) %>%
      dplyr::select(feature1, feature2, pearson_corr)

    ## Name nodes
    graph_table <- plot_data %>%
      tidygraph::as_tbl_graph() %>%
      dplyr::mutate(type = dplyr::case_when(name %in% input$content ~ "Query",
                                            name %in% data_universal_achilles_long$gene ~ "Gene",
                                            !(name %in% data_universal_achilles_long$gene) ~ "Pathway")
      )

    ## Plot network
    plot_complete <- ggraph::ggraph(graph_table, layout = "fr") +
      ggraph::geom_edge_link(ggplot2::aes(edge_alpha = abs(pearson_corr)),
                             edge_width = 0.5, show.legend = FALSE) +
      ggraph::geom_node_point(ggplot2::aes(fill = type), color = "black",
                              pch = 21, size = 3, alpha = 0.85) +
      {if(!is.null(highlight))ggraph::geom_node_point(ggplot2::aes(filter = name %in% highlight),
                                                      fill = "red",
                                                      pch = 21,
                                                      size = 3,
                                                      alpha = 0.85)} +
      {if(!is.null(highlight) & show_labels)ggraph::geom_node_label(ggplot2::aes(label = name, filter = name %in% highlight),
                                                                    repel = TRUE,
                                                                    size = fontsize,
                                                                    fill = ggplot2::alpha(c("white"), 0.6),
                                                                    label.size = NA,
                                                                    fontface = "bold",
                                                                    family = "Roboto Slab")} +
      ggraph::theme_graph(foreground = "white", fg_text_colour = "white") +
      ggplot2::theme(legend.position = "top",
                     legend.title = ggplot2::element_blank(),
                     text = ggplot2::element_text(family = "Roboto Slab", size = 16)) +
      scale_fill_ddh_d() +
      ggraph::scale_edge_color_continuous(low = "black", high = "black")

    return(plot_complete)
  }

  #error handling
  tryCatch(make_gene_pathways_components_network_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})

}

## EXPvDEP PLOT --------------------------------------------------------
#' Gene Dependency versus Expression
#'
#' Each point shows the dependency value compared to the expression value for gene within a given cell line. Gray area indicates dependency values that are between -1 and 1.
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a expdep plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_expdep(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_expdep(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_expdep(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_expdep <- function(data_universal_expression_long = universal_expression_long,
                        data_universal_achilles_long = universal_achilles_long,
                        data_cell_expression_names = cell_expression_names,
                        plot_se = TRUE,
                        input = list(),
                        card = FALSE) {
  make_expdep_raw <- function() {

    if (input$type == "gene") {
      exp_data <-
        data_universal_expression_long %>% #plot setup
        dplyr::select(X1, gene, gene_expression) %>%
        dplyr::filter(gene %in% input$content)

      dep_data <-
        data_universal_achilles_long %>% #plot setup
        dplyr::filter(gene %in% input$content)

      combined_data <-
        exp_data %>%
        dplyr::inner_join(dep_data, by = c("X1", "gene")) %>%
        dplyr::filter(!is.na(dep_score),
                      !is.na(gene_expression)) %>%
        dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
        dplyr::mutate(med = median(dep_score, na.rm = TRUE)) %>%
        dplyr::left_join(data_cell_expression_names, by = "X1") %>%
        dplyr::select(gene, gene_expression, dep_score, med, cell_line, lineage)
    } else if (input$type == "cell") {
      exp_data <-
        data_universal_expression_long %>% #plot setup
        dplyr::left_join(data_cell_expression_names, by = "X1") %>%
        dplyr::select(gene, cell_line, gene_expression) %>%
        dplyr::filter(cell_line %in% input$content)

      dep_data <-
        data_universal_achilles_long %>% #plot setup
        dplyr::left_join(data_cell_expression_names, by = "X1") %>%
        dplyr::filter(cell_line %in% input$content)

      combined_data <-
        exp_data %>%
        dplyr::inner_join(dep_data, by = c("cell_line", "gene")) %>%
        dplyr::filter(!is.na(dep_score),
                      !is.na(gene_expression)) %>%
        dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
        dplyr::mutate(med = median(dep_score, na.rm = TRUE)) %>%
        dplyr::select(gene, gene_expression, dep_score, med, cell_line, lineage)
    }

    plot_complete <-
      combined_data %>%
      ggplot2::ggplot(ggplot2::aes(dep_score, gene_expression)) +
      ## gray background
      ## dot plot
      {if(input$type == "gene")ggplot2::geom_point(ggplot2::aes(color = forcats::fct_reorder(gene, med),
                                                                fill = forcats::fct_reorder(gene, med)),
                                                   size = 2, stroke = .1, alpha = 0.4)} +
      {if(input$type == "cell")ggplot2::geom_point(ggplot2::aes(color = forcats::fct_reorder(cell_line, med),
                                                                fill = forcats::fct_reorder(cell_line, med)),
                                                   size = 2, stroke = .1, alpha = 0.4)} +
      # smooth line
      {if(input$type == "gene")ggplot2::geom_smooth(ggplot2::aes(color = forcats::fct_reorder(gene, med),
                                                                 fill = forcats::fct_reorder(gene, med)),
                                                    method = "lm",
                                                    se = plot_se)} +
      {if(input$type == "cell")ggplot2::geom_smooth(ggplot2::aes(color = forcats::fct_reorder(cell_line, med),
                                                                 fill = forcats::fct_reorder(cell_line, med)),
                                                    method = "lm",
                                                    se = plot_se)} +
      # R coefs
      {if(card == FALSE & input$type == "gene")ggpubr::stat_cor(ggplot2::aes(color = forcats::fct_reorder(gene, med)),
                                                                digits = 3)} +
      {if(card == FALSE & input$type == "cell")ggpubr::stat_cor(ggplot2::aes(color = forcats::fct_reorder(cell_line, med)),
                                                                digits = 3)} +
      ## scales + legends
      scale_color_ddh_d(palette = input$type) +
      scale_fill_ddh_d(palette = input$type) +
      ggplot2::guides(
        color = ggplot2::guide_legend(reverse = TRUE, override.aes = list(size = 5)),
        fill = ggplot2::guide_legend(reverse = TRUE, override.aes = list(size = 5))
      ) +
      ## titles
      ggplot2::labs(
        x = "Dependency",
        y = "Expression",
        color = ifelse(input$type == "gene", "Query Gene", "Query Cell Line"),
        fill = ifelse(input$type == "gene", "Query Gene", "Query Cell Line")
      ) +
      ## theme changes
      theme_ddh() +
      ggplot2::theme(
        text = ggplot2::element_text(family = "Nunito Sans"),
        axis.text = ggplot2::element_text(family = "Roboto Slab")
      ) +
      NULL

    if(length(input$content) == 1){
      plot_complete  <-
        plot_complete +
        ggplot2::guides(color = "none",
                        fill = "none") +
        ggplot2::labs(x = paste0(input$content, " Dependency"),
                      y = paste0(input$content, " Expression"))
    }

    if(card == TRUE) {
      plot_complete <-
        plot_complete +
        ggplot2::labs(x = "", y = "") +
        ggplot2::theme(legend.position='none') +
        NULL
    }
    return(plot_complete)
  }
  #error handling
  tryCatch(make_expdep_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

#CELL --------------------------------------------------------------------
#' Cell Images
#'
#' Cells image shown at low (top) and high (bottom) growth density.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_cell_image(input = list(content = "HEPG2"))
make_cell_image <- function(input = list(),
                            card = FALSE) {
  make_cell_image_raw <- function(){
    if(length(input$content > 1)){
      cell_name <- sample(input$content, 1) #input$content[[1]]
    } else {
      cell_name <- input$content
    }

    #image type
    if(card == TRUE){image_type = "card"} else {image_type = "plot"}

    #check if exists
    url <- glue::glue("https://{Sys.getenv('AWS_CELLIMAGES_BUCKET_ID')}.s3.amazonaws.com/{cell_name}_cell_image_{image_type}.jpeg")
    status <- httr::GET(url) %>% httr::status_code()

    if(status == 200){
      return(url)
    } else {
      num <- sample(1:5, 1)
      return(glue::glue("https://{Sys.getenv('AWS_CELLIMAGES_BUCKET_ID')}.s3.amazonaws.com/error_cell_image{num}.jpg"))
    }
  }
  #error handling
  tryCatch(make_cell_image_raw(),
           error = function(e){
             message(e)
             num <- sample(1:5, 1)
             return(glue::glue("https://{Sys.getenv('AWS_CELLIMAGES_BUCKET_ID')}.s3.amazonaws.com/error_cell_image{num}.jpg"))
           })
}

## CO-ESSENTIALITY CELL LINE PLOT --------------------------------------------------------
#' Cell Similarity Plot
#'
#' Each point shows the ranked linear model coefficient estimate value ordered from high to low for each query.
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a cell similarity plot If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_cell_similarity(input = list(type = "cell", query = "HEPG2", content = "HEPG2"))
#' make_cell_similarity(input = list(type = "cell", query = "HEPG2", content = "HEPG2"), similarity = "expression")
#' \dontrun{
#' make_cell_similarity(input = list(type = "cell", query = "HEPG2", content = "HEPG2"))
#' }
make_cell_similarity <- function(data_cell_dependency_sim = cell_dependency_sim,
                                 data_cell_expression_sim = cell_expression_sim,
                                 similarity = "dependency",
                                 input = list(),
                                 card = FALSE,
                                 scale = NULL) { #no card option, but need this to prevent error
  make_similarity_raw <- function() {

    if(similarity == "dependency") {
      cell_sims <- data_cell_dependency_sim
    } else if(similarity == "expression") {
      cell_sims <- data_cell_expression_sim
    }

    cell_sim_table <-
      cell_sims %>%
      dplyr::filter(cell1_name %in% input$content | cell2_name %in% input$content)

    # Swap cols
    cell_sim_table[cell_sim_table$cell2_name %in% input$content, c("cell1_name", "cell2_name")] <-
      cell_sim_table[cell_sim_table$cell2_name %in% input$content, c("cell2_name", "cell1_name")]

    content_var <- glue::glue_collapse(input$content, sep = ", ")

    plot_data <-
      cell_sim_table %>%
      dplyr::group_by(cell1_name) %>%
      dplyr::arrange(dplyr::desc(coef)) %>%
      dplyr::mutate(
        rank = 1:dplyr::n(),
        med = median(coef, na.rm = TRUE)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::rename(name = cell1_name) #for plot name


    if(!is.null(scale) & !card){
      plot_data <-
        plot_data %>%
        dplyr::slice_sample(prop = scale)
    }

    if(card) {
      if(is.null(scale)) {
        scale <- 0.3
      }
      plot_data <- plot_data %>%
        dplyr::group_by(name) %>%
        dplyr::sample_n(scale*dplyr::n()) %>% # scale is also used here
        dplyr::ungroup()
    }

    plot_complete <-
      plot_data %>%
      ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = 0, color = "gray80", linetype = "dashed") +
      ggplot2::annotate("text", x = Inf, y = 0.005, label = "Cell Rank",
                        color = "gray40", hjust = 1.15 ,vjust = 0) +
      ## dot plot
      ggplot2::geom_point(ggplot2::aes(x = rank,
                                       y = coef,
                                       text = glue::glue('Cell: {name}'),
                                       color = forcats::fct_reorder(name, med),
                                       fill = forcats::fct_reorder(name, med)
      ),
      size = 1.1, stroke = .1, alpha = 0.4) +
      ## scales + legends
      ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = 0.02), na.translate = FALSE) +
      scale_color_ddh_d(palette = input$type) +
      scale_fill_ddh_d(palette = input$type) +
      ggplot2::guides(
        color = ggplot2::guide_legend(reverse = TRUE, override.aes = list(size = 5)),
        fill = ggplot2::guide_legend(reverse = TRUE, override.aes = list(size = 5))
      ) +
      ## labels
      ggplot2::labs(
        x = NULL,
        y = glue::glue("Cell coefficient estimates with {content_var}"),
        color = "Query Cell",
        fill = "Query Cell"
      ) +
      ## theme changes
      theme_ddh() +
      ggplot2::theme(
        text = ggplot2::element_text(family = "Nunito Sans"),
        axis.text = ggplot2::element_text(family = "Roboto Slab"),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank()
      ) +
      NULL

    ##change axis in case of more than 1 gene selected
    if(length(input$content) == 1){ #fix me
      plot_complete  <-
        plot_complete +
        ggplot2::guides(color = "none",
                        fill = "none")
    }

    return(plot_complete)
  }
  #error handling
  tryCatch(make_similarity_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

# FUNCTIONAL PLOT ------------------------------------------
#' Differential Pathway Expression Plot
#'
#' The colored vertical bars indicate the pathway median expression for the queried cell line/s while the background grey points indicate the pathway median expression of all the other cell lines. If the query includes only one cell line, the difference between the median pathway expression of that cell line and the median pathway expression of all the other cell lines will be computed and the gene_pathways with higher differences will appear first in the plot. Otherwise, the biggest difference between pathway medians of the queried cell lines will be used to rank the gene_pathways in the plot. Those gene_pathways with higher differences will appear first in the plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_functional_cell(input = list(type = "cell", query = "HEPG2", content = "HEPG2"))
#' make_functional_cell(input = list(type = "cell", query = c("HEPG2", "HEL"), content = c("HEPG2", "HEL")))
#' make_functional_cell(input = list(type = "cell", query = "HEPG2", content = "HEPG2"), card = TRUE)
#' \dontrun{
#' make_functional_cell(input = list(type = "cell", query = "HEPG2", content = "HEPG2"))
#' }
make_functional_cell <- function(data_gene_pathways = gene_pathways,
                                 data_universal_expression_long = universal_expression_long,
                                 data_cell_expression_meta = cell_expression_meta,
                                 input = list(),
                                 num_pathways = 10,
                                 nwords = 5,
                                 num_genes = 2,
                                 remove_equivalent_pathways = FALSE,
                                 card = FALSE) {
  make_functional_cell_raw <- function() {

    plot_data <-
      data_gene_pathways %>%
      tidyr::unnest(data) %>%
      dplyr::left_join(data_universal_expression_long, by = "gene") %>%
      dplyr::select(pathway, go, gene, gene_expression, X1) %>%
      tidyr::drop_na() %>%
      dplyr::left_join(data_cell_expression_meta %>%
                         dplyr::select(X1, cell_line), by = "X1") %>%
      dplyr::select(-X1) %>%
      dplyr::mutate(pathway_short = ifelse(stringr::str_count(pathway) >= nwords,
                                           paste0(gsub(paste0("^((\\w+\\W+){", nwords, "}\\w+).*$"), "\\1", pathway), " ..."),
                                           pathway)
      )

    med_cell <-
      plot_data %>%
      dplyr::filter(cell_line %in% input$content) %>%
      dplyr::group_by(go, cell_line) %>%
      dplyr::mutate(genes_num = dplyr::n()) %>%
      dplyr::filter(genes_num >= num_genes) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(cell_line, go) %>%
      dplyr::summarise(med = median(gene_expression)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(med != 0) %>% # decide if filter zero expressions
      dplyr::left_join(plot_data %>%
                         dplyr::select(go, pathway_short) %>%
                         dplyr::filter(!duplicated(go)),
                       by = "go")

    if(remove_equivalent_pathways) {
      med_cell <-
        med_cell %>%
        dplyr::filter(!duplicated(med))
    }

    if(length(input$content) == 1) {
      diff_cell <-
        plot_data %>%
        dplyr::filter(!cell_line %in% input$content) %>%
        dplyr::group_by(go) %>%
        dplyr::mutate(genes_num = dplyr::n()) %>%
        dplyr::filter(genes_num > num_genes) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(go) %>%
        dplyr::summarise(med = median(gene_expression)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(med != 0) %>% # decide if filter zero expressions
        dplyr::bind_rows(med_cell) %>%
        dplyr::group_by(go) %>%
        dplyr::summarise(diff = abs(diff(med))) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(-diff) %>%
        dplyr::slice(1:num_pathways) %>%
        dplyr::left_join(plot_data %>%
                           dplyr::select(go, pathway_short) %>%
                           dplyr::filter(!duplicated(go)),
                         by = "go")
    } else {
      diff_cell <-
        med_cell %>%
        dplyr::group_by(go) %>%
        dplyr::summarise(diff = abs(diff(med))) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(-diff) %>%
        dplyr::slice(1:num_pathways) %>%
        dplyr::left_join(plot_data %>%
                           dplyr::select(go, pathway_short) %>%
                           dplyr::filter(!duplicated(go)),
                         by = "go")
    }

    med_cell <-
      med_cell %>%
      dplyr::filter(go %in% diff_cell$go) %>%
      dplyr::left_join(diff_cell %>%
                         dplyr::select(-pathway_short) %>%
                         dplyr::filter(!duplicated(go)),
                       by = "go")

    background <-
      plot_data %>%
      dplyr::filter(go %in% diff_cell$go) %>%
      dplyr::filter(!cell_line %in% input$content) %>%
      dplyr::group_by(cell_line, go) %>%
      dplyr::summarise(med = median(gene_expression)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(med != 0) %>% # decide if filter zero expressions
      dplyr::left_join(diff_cell, by = "go")

    plot_complete <-
      ggplot2::ggplot() +
      ggplot2::geom_jitter(data = background, ggplot2::aes(reorder(pathway_short, diff), med), color = "gray69", alpha = 0.5, width = 0.25) +
      ggplot2::geom_crossbar(data = med_cell,
                             ggplot2::aes(reorder(pathway_short, diff), med, ymin = med, ymax = med, color = cell_line)) +
      ggplot2::geom_point(data = med_cell, ggplot2::aes(reorder(pathway_short, diff), med, color = cell_line), size = 3) +
      ggplot2::labs(x = NULL,
                    y = "Gene Expression",
                    color = "Query Cell") +
      ggplot2::coord_flip() +
      scale_color_ddh_d(palette = input$type) +
      ggplot2::guides(
        color = ggplot2::guide_legend(reverse = TRUE, override.aes = list(size = 5))
      ) +
      ## theme changes
      theme_ddh() +
      ggplot2::theme(
        text = ggplot2::element_text(family = "Nunito Sans"),
        axis.text = ggplot2::element_text(family = "Roboto Slab"),
        legend.title = ggplot2::element_blank()
      ) +
      NULL

    if(length(input$content) == 1){
      plot_complete  <-
        plot_complete +
        ggplot2::ylab(paste0(stringr::str_c(input$content, collapse = ", "), " Gene Expression")) +
        ggplot2::theme(
          legend.position = "none"
        )
    } else {
      plot_complete <- plot_complete +
        ggplot2::theme(legend.title = ggplot2::element_blank())
    }

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_blank(),
                       legend.position = "none"
        )
    }

    return(plot_complete)
  }

  #error handling
  tryCatch(make_functional_cell_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

# CELL METADATA PLOT ------------------------------------------
#' Lineage Similarity Plot
#'
#' Similar and dissimilar lineages (and sublineages) associated with the queried cell line/s.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_metadata_cell(input = list(type = "cell", content = c("HEPG2")))
#' make_metadata_cell(input = list(type = "cell", content = c("HEPG2")), cell_line_similarity = "expression")
#' make_metadata_cell(input = list(type = "cell", content = c("HEPG2")), metadata = "sublineage")
#' make_metadata_cell(input = list(type = "cell", content = c("HEPG2", "HEL")))
#' make_metadata_cell(input = list(type = "cell", content = c("HEPG2")), card = TRUE)
#' \dontrun{
#' make_metadata_cell(input = list(type = "cell", content = c("HEPG2")))
#' }
make_metadata_cell <- function(data_cell_dependency_sim = cell_dependency_sim,
                               data_cell_expression_sim = cell_expression_sim,
                               input = list(),
                               cell_line_similarity = "dependency",
                               metadata = "lineage",
                               bonferroni_cutoff = 0.05,
                               card = FALSE) {
  make_metadata_cell_raw <- function() {

    plot_data <- make_cell_sim_table(data_cell_dependency_sim,
                                     data_cell_expression_sim,
                                     similarity = cell_line_similarity,
                                     bonferroni_cutoff = 1.1, # to include bonf == 1
                                     input = input) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(group = dplyr::case_when(bonferroni > bonferroni_cutoff ~ "None",
                                             bonferroni < bonferroni_cutoff & coef > 0 ~ "Similar",
                                             bonferroni < bonferroni_cutoff & coef < 0 ~ "Dissimilar"),
                    group = forcats::as_factor(group)
      ) %>%
      dplyr::as_tibble()

    if(metadata == "lineage") {
      siglin <-
        plot_data %>%
        dplyr::filter(bonferroni < bonferroni_cutoff) %>%
        dplyr::pull(lineage) %>%
        unique()

      plot_data <-
        plot_data %>%
        dplyr::filter(lineage %in% siglin)

      if(card & nrow(plot_data) > 6) {
        most_associated <-
          plot_data %>%
          dplyr::arrange(pval) %>%
          dplyr::filter(!duplicated(lineage)) %>%
          dplyr::slice(1:6) %>%
          dplyr::pull(lineage)

        plot_data <-
          plot_data %>%
          dplyr::filter(., lineage %in% most_associated)
      }

      plot_complete <-
        plot_data %>%
        dplyr::filter(!is.na(lineage)) %>%
        dplyr::group_by(lineage, group) %>%
        dplyr::count() %>%
        dplyr::ungroup() %>%
        ggplot2::ggplot(ggplot2::aes(x = n, y = lineage, fill = group, group = group))
    }
    else if (metadata == "sublineage") {
      sigsublin <-
        plot_data %>%
        dplyr::filter(bonferroni < bonferroni_cutoff) %>%
        dplyr::pull(lineage_subtype) %>%
        unique()

      plot_data <-
        plot_data %>%
        dplyr::filter(lineage_subtype %in% sigsublin)

      plot_complete <-
        plot_data %>%
        dplyr::filter(!is.na(lineage_subtype)) %>%
        dplyr::group_by(lineage_subtype, group) %>%
        dplyr::count() %>%
        dplyr::ungroup() %>%
        ggplot2::ggplot(ggplot2::aes(x = n, y = lineage_subtype, fill = group, group = group))
    }

    plot_complete <-
      plot_complete +
      ggplot2::geom_col() +
      ggplot2::labs(x = "# Cell Lines",
                    y = NULL) +
      # scale_x_continuous(labels = scales::percent_format()) +
      scale_fill_ddh_d(palette = input$type) +
      ggplot2::guides(
        color = ggplot2::guide_legend(reverse = TRUE, override.aes = list(size = 5))
      ) +
      ## theme changes
      theme_ddh() +
      ggplot2::theme(
        text = ggplot2::element_text(family = "Nunito Sans"),
        axis.text = ggplot2::element_text(family = "Roboto Slab"),
        legend.title = ggplot2::element_blank(),
        legend.position = "top"
      ) +
      NULL

    if(card) {
      plot_complete <-
        plot_complete +
        ggplot2::scale_y_discrete(labels = scales::label_wrap(10)) + #to prevent long lines and squished plots
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       plot.title = ggplot2::element_blank(),
                       legend.position = "none"
        )
    }
    return(plot_complete)
  }

  #error handling
  tryCatch(make_metadata_cell_raw(),
           error = function(e){
             message(e)
             make_bomb_plot()})
}

#COMPOUND --------------------------------------------------------------------
#' Molecule Structure Plot
#'
#' \code{make_molecule_structure} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a molecule structure plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a molecule structure plot. If an error is thrown, then will return a bomb plot.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_molecule_structure(input = list(content = "aspirin"))
#' make_molecule_structure(input = list(content = "ibuprofen"))
#' make_molecule_structure(input = list(content = 2244))
#' \dontrun{
#' make_molecule_structure(input = list(content = "aspirin"))
#' }
make_molecule_structure <- function(input = list(),
                                    # file_name = NULL,
                                    card = FALSE) {
  make_molecule_structure_raw <- function(){
    #Customize the material with toon shading
    shiny_toon_material = rayvertex::material_list(type="toon_phong",
                                                   toon_levels=,
                                                   toon_outline_width=0.1)

    molecule <-
      raymolecule::get_molecule(input$content) %>%
      raymolecule::generate_full_scene(pathtrace=FALSE,
                                       material_vertex = shiny_toon_material)
    if(card == TRUE){
      model_width = 1080
      model_height = 1080
    } else {
      model_width = 2160
      model_height = 2160
    }
    # if(!is.null(file_name)){
    #   molecule %>%
    #     raymolecule::render_model(width=model_width,
    #                               height=model_height,
    #                               fov=10,
    #                               background="white",
    #                               filename = file_name)
    #   return(message(glue::glue('{input$content} saved to {file_name}')))
    # } else {
      molecule %>%
        raymolecule::render_model(width=model_width,
                                  height=model_height,
                                  background="white")
    # }
  }
  #error handling
  tryCatch(make_molecule_structure_raw(),
           error = function(e){
             message(e)
           })
}
