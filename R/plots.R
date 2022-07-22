
## BARCODE PLOT --------------------------------------------------------------------
make_barcode <- function(input = list()) {
  make_barcode_raw <- function() {
    #if multiple, then pull single "rep" image; consider pulling >1 and using patchwork, eg.
    if(length(input$content > 1)){
      fav_gene <- sample(input$content, 1) #input$content[[1]]
    } else {
      fav_gene <- input$content
    }

    #fetch image from dir, HARDCODED TO DATA DIR, NOT TEST DATA DIR
    file_name <- glue::glue('{fav_gene}_barcode_card.jpeg')
    image_path <- here::here("data", "images", "gene", fav_gene, file_name)

    # barcode_image <-
    #   magick::image_read(path = image_path)

    #error catching
    #currently tryCatch will just make empty image

    #place image on ggplot
    # plot_complete <-
    #   ggplot() +
    #   ggpubr::background_image(barcode_image) +
    #   theme_void()

    return(image_path)
  }
  #error handling
  tryCatch(make_barcode_raw(),
           error = function(x){make_bomb_plot()})
}

## IDEOGRAM PLOT --------------------------------------------------------
#' Ideogram Plot
#'
#' \code{make_ideogram} returns an image of a chromosome with a gene position annotated.
#'
#' This is a plotting function that takes a gene name and returns
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns an ideogram plot. If an error is thrown, then will return a bomb plot

#' @export
#' @examples
#' sum(1:10)
#' make_ideogram(input = list(type = "gene", query = "ROCK1", content = "ROCK1"))
#' make_ideogram(input = list(type = "gene", query = "ROCK1", content = "ROCK1"), card = TRUE)
#' \dontrun{
#' make_ideogram(input = list(type = "gene", content = "ROCK1"))
#' }
make_ideogram <- function(location_data = gene_location,
                          chromosome_data = chromosome,
                          input = list(),
                          card = FALSE) {
  make_ideogram_raw <- function() {
    #set baseline chromosome info
    chromosome_list <- pull(chromosome, id)
    chromosome_levels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
    chromosome_list <- fct_relevel(chromosome_list, chromosome_levels)

    #adjust by n so it bumps genes/bands off chromosome ends
    n <- 2000000 #chr1 is 250M

    #get some pq arms for the background images
    pq <-
      chromosome_data %>%
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
      chromosome_data %>%
      dplyr::select(chromosome_name = id, basepairs) %>%
      dplyr::mutate(basepairs = basepairs + n) #correct for n adjustment here

    #make some bands
    band_boundaries <-
      location_data %>%
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
      dplyr::mutate(alternate = row_number(chromosome_name) %% 2)

    #adding if (is.null(gene_symbol)), gene_symbol=NULL allows me to generate a plot, even without gene_symbol data, #no data filters, #drop geom_point here
    #get location data for each gene query
    gene_loci <-
      location_data %>%
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
      filter(chromosome_name %in% chromosome_loci) %>%
      slice_max(basepairs) %>%
      pull(basepairs)

    ideogram_plot <-
      ggplot() +
      #background line fixes height
      geom_segment(data = pq %>% filter(chromosome_name %in% chromosome_loci),
                   aes(x = chromosome_name, xend = chromosome_name, y = 0, yend = 260000000, alpha = 1),
                   color = "white", size = 1, lineend = "butt") +
      #background for black line
      geom_segment(data = pq %>% filter(chromosome_name %in% chromosome_loci),
                   aes(x = chromosome_name, xend = chromosome_name, y = y_start, yend = y_end, alpha = 0.5),
                   color = "black", size = 7, lineend = "round") +
      #chromosome
      geom_segment(data = pq %>% filter(chromosome_name %in% chromosome_loci),
                   aes(x = chromosome_name, xend = chromosome_name, y = y_start, yend = y_end),
                   color = "gray95", size = 6, lineend = "round") +
      #centromere
      geom_point(data = pq %>% filter(chromosome_name %in% chromosome_loci,
                                      arm == "q"),
                 aes(x = chromosome_name, y = y_end + n), #calculated 1/2 of distance between
                 color = "gray90", size = 5.5) +
      #pq bands
      geom_segment(data = band_boundaries %>% filter(alternate == 1, chromosome_name %in% chromosome_loci),
                   aes(x = chromosome_name, xend = chromosome_name, y = min, yend = max),
                   size = 6.5, color = "black", alpha = 0.5) +
      #gene points + labels
      geom_point(data = gene_loci, aes(x = chromosome_name, y = start),  size = 6, color = ddh_pal_d(palette = "gene")(1), alpha = 1) +
      ggrepel::geom_text_repel(data = gene_loci, aes(x = chromosome_name, y = start, label = approved_symbol), nudge_x = .2, min.segment.length = 1, family = "Chivo") +
      scale_fill_manual(values = c("#FFFFFF", "#FFFFFF")) +
      labs(y = NULL) +
      theme_ddh(base_size = 16) +
      theme_void() +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 12)) +
      #clip height
      coord_cartesian(ylim=c(0, clip_height)) +
      NULL

    if(card == TRUE){
      ideogram_plot <-
        ideogram_plot +
        labs(x = "") + #, title = "Gene Information", caption = "more ...") +
        NULL
    }
    return(ideogram_plot)
  }
  #error handling
  tryCatch(make_ideogram_raw(),
           error = function(x){make_bomb_plot()})
}

#figure legend
plot_ideogram_title <- "Gene Location."
plot_ideogram_legend <- paste0("Each point shows the location of the query gene(s) on human chromosomes.")

## SIZE PLOT --------------------------------------------------------
#' Protein Size Plot
#'
#' \code{make_proteinsize} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a proteinsize plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a proteinsize plot. If an error is thrown, then will return a bomb plot.
#'
#' @export
#' @examples
#' make_proteinsize(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_proteinsize(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_proteinsize(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_proteinsize <- function(protein_data = proteins,
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
      proteins %>%
      dplyr::filter(gene_name %in% input$content) %>%
      dplyr::arrange(mass) %>%
      dplyr::pull(gene_name)

    moreThanThree <- length(gene_symbol) > 3
    moreThanFour<- length(gene_symbol) > 4

    last_gene <- gene_symbol[length(gene_symbol)]
    last_mass <-
      protein_data %>%
      dplyr::filter(gene_name == last_gene) %>%
      dplyr::pull(mass)

    make_mass_strip <- function(data = proteins,
                                var,
                                card_var = card,
                                max_mass = last_mass) {

      selected <-
        protein_data %>%
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
        proteins %>%
        ggplot(aes(x = mass)) +
        ## draw distribution as colored strip
        ggdist::stat_interval(
          aes(y = 1),
          .width = w, size = 10
        ) +
        ## line locator
        geom_linerange(
          data = selected, aes(ymin = .9, ymax = 1.03),
          color = "white", size = .8
        ) +
        ## add triangular locator symbol
        geom_point(
          data = selected, aes(y = triangle_y),
          shape = 6, size = triangle_size, stroke = 1
        ) +
        ## add gene name label
        {if(!card_var)geom_text( #| !moreThanFour
          data = selected, aes(y = gene_y, label = gene_name),
          family = "Chivo", size = 5.2, vjust = 0
        )} +
        ## add mass label
        {if(!card_var)geom_text( # | !moreThanThree
          data = selected, aes(y = mass_y, label = paste(mass, "kDa")),
          family = "Chivo", size = 4.3, vjust = 0
        )} +
        scale_x_continuous(limits = c(0, max_mass*2),
                           labels = function(x) paste(x, "kDa")) +
        coord_cartesian(ylim = c(.95, 1.2)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_color_manual(values = colors, guide = "none") +
        theme_ddh() +
        theme_void() +
        NULL


      if (var == last_gene) {
        base_plot <-
          base_plot +
          labs(x = "Protein Size") +
          theme(axis.text.x = element_text(family = "Roboto Slab", size = 14),
                axis.title.x = element_text(family = "Nunito Sans", size = 18,
                                            margin = margin(t = 12, b = 12)))
      }

      return(base_plot)
    }

    glist <- lapply(gene_symbol, make_mass_strip, data = proteins)
    plot_complete <- patchwork::wrap_plots(glist, nrow = length(glist))

    if(card == TRUE){
      plot_complete <- make_mass_strip(var = gene_symbol)

      mass <- #get mass to center clipping in next step
        proteins %>%
        filter(gene_name %in% input$content) %>%
        pull(mass) %>%
        median(na.rm = TRUE)

      plot_complete <-
        plot_complete +
        coord_cartesian(ylim = c(.9, 1.1), xlim = c(mass - mass*.3, mass + mass*.3)) +
        # plot_annotation(
        #   title = "Size Information",
        #   caption = "more ...") +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank()
        )
    }

    return(plot_complete)
  }
  #error handling
  tryCatch(make_proteinsize_raw(),
           error = function(x){make_bomb_plot()})
}

#test
# make_proteinsize(input = list(content = "ROCK1"))
# make_proteinsize(input = list(content = "ROCK1"), card = TRUE)
# make_proteinsize(input = list(content = c("ROCK1", "ROCK2")))
# make_proteinsize(input = list(content = c("ROCK1", "ROCK2")), card = TRUE)
# make_proteinsize(input = list(content = c("RDX", "ROCK2", "DTX3L", "MSN", "SORL1", "EZR")))
# make_proteinsize(input = list(content = c("RDX", "ROCK2", "DTX3L", "MSN", "SORL1", "EZR")), card = TRUE)

#figure legend
plot_proteinsize_title <- "Protein Size."
plot_proteinsize_legend <- paste0("Mass compared to all protein masses. The colored strip visualizes the distribution of protein sizes. Each colored box is thus representing a decile of the full data. The triangle indicates where exactly the queried genes fall on this gradient of protein sizes.")

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
#' @export
#' @examples
#' make_sequence(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_sequence(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_sequence(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_sequence <- function(sequence_data = proteins,
                          input = list(),
                          card = FALSE) {
  make_sequence_raw <- function() {
    #if multiple, then pull single "rep" image; consider pulling >1 and using patchwork, eg.
    # if(length(input$content > 1)){
    #   fav_gene <- input$content[[1]] #sample(input$content, 1)
    # } else {
    #   fav_gene <- input$content
    # }

    sequence_string <-
      sequence_data %>%
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
      map2(.x = sequence_vec,
           .y = sequence_num,
           .f = get_subsequence)

    sequence_font <- "Roboto Mono" #"Roboto Slab"
    sequence_size <- 10

    plot_complete <-
      ggplot() +
      ## annotate with text
      annotate("text", x = 0, y = 10, label = sequence_complete[1], family = sequence_font, size = sequence_size) +
      annotate("text", x = 0, y = 9, label = sequence_complete[2], family = sequence_font, size = sequence_size) +
      annotate("text", x = 0, y = 8, label = sequence_complete[3], family = sequence_font, size = sequence_size) +
      annotate("text", x = 0, y = 7, label = sequence_complete[4], family = sequence_font, size = sequence_size) +
      annotate("text", x = 0, y = 6, label = sequence_complete[5], family = sequence_font, size = sequence_size) +
      annotate("text", x = 0, y = 5, label = sequence_complete[6], family = sequence_font, size = sequence_size) +
      annotate("text", x = 0, y = 4, label = sequence_complete[7], family = sequence_font, size = sequence_size) +
      annotate("text", x = 0, y = 3, label = sequence_complete[8], family = sequence_font, size = sequence_size) +
      annotate("text", x = 0, y = 2, label = sequence_complete[9], family = sequence_font, size = sequence_size) +
      annotate("text", x = 0, y = 1, label = sequence_complete[10], family = sequence_font, size = sequence_size) +
      coord_cartesian(ylim = c(0.5, 10.5), clip = 'on') +
      theme_void() +
      NULL

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        labs(x = "") +#, title = "Sequence Information", caption = "more ...") +
        NULL
    }

    return(plot_complete)
  }
  #error handling
  tryCatch(make_sequence_raw(),
           error = function(x){make_bomb_plot()})
}

#make_sequence(input = list(content = "ROCK2"))
#make_sequence(input = list(content = c("ROCK1", "ROCK2")))
#make_sequence(input = list(content = "ROCK2"), card = TRUE)

## PROTEIN DOMAIN PLOT --------------------------------------------------------
#' Protein Domain Plot
#'
#' \code{make_protein_domain_plot} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a protein domain plot plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a protein domain plot. If an error is thrown, then will return a bomb plot.
#'
#' @export
#' @examples
#' make_protein_domain_plot(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_protein_domain_plot(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_protein_domain_plot(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_protein_domain_plot <- function(input = list(),
                                     domain_data = protein_domains,
                                     dom_var = NULL,
                                     ptm_var = NULL) {

  make_protein_domain_plot_raw <- function() {

    gene_symbol <- domain_data %>%
      filter(gene_name %in% input$content) %>%
      pull(gene_name) %>%
      unique()

    lengths_data <- domain_data %>%
      filter(gene_name %in% input$content) %>%
      filter(!duplicated(gene_name))

    plot_data <- domain_data %>%
      dplyr::filter(gene_name %in% input$content) %>%
      mutate(order = 1)

    prots_dr <- plot_data %>%
      filter(type %in% c("DOMAIN", "REGION")) %>%
      drop_na(description)

    prots_ptm <- plot_data %>%
      filter(category == "PTM") %>%
      drop_na(description)

    # Plot
    base_plot <-
      ggplot() +
      geom_segment(data = plot_data,
                   aes(x = 1,
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
        geom_rect(data = prots_dr %>%
                    filter(description %in% dom_var),
                  aes(xmin = begin,
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
        geom_point(data = prots_ptm %>%
                     filter(description %in% ptm_var),
                   aes(x = begin,
                       y = order + 3,
                       shape = description),
                   size = 8,
                   alpha = 1,
                   color = "black") +
        geom_segment(data = prots_ptm %>%
                       filter(description %in% ptm_var),
                     aes(x = begin,
                         xend = begin,
                         y = order,
                         yend = order + 3),
                     size = 0.5,
                     linetype = 2,
                     color = "black") +
        scale_y_continuous(limits=c(0,4.1)) + #should be a touch larger than yend
        NULL
    }

    wrapped_labels <- function(.seq_len) {
      fun <- function(labels) {
        labels <- label_value(labels, multi_line = TRUE)
        lapply(labels, function(x) {
          x <- paste0(x, " (", .seq_len, " amino acids)")
          vapply(x, paste, character(1), collapse = "\n")
        })
      }
      structure(fun, class = "labeller")
    }

    base_plot <-
      base_plot +
      labs(x = NULL,
           y = NULL) +
      facet_wrap(~ gene_name, ncol = 1,
                 labeller = wrapped_labels(.seq_len = lengths_data$seq_len)
      ) +
      theme_ddh(base_size = 16) +
      theme_void() +
      theme(
        legend.title = element_blank(),
        strip.text = element_text(family = "Chivo", hjust = 0.5, size = 15),
        plot.title = element_text(family = "Chivo", hjust = 0.5, size = 15),
        text = element_text(family = "Nunito Sans"),
        legend.text = element_text(size = 15),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      ) +
      NULL

    plot_complete <-
      base_plot +
      theme(axis.text.x = element_text(family = "Roboto Slab", size = 14),
            axis.title.x = element_text(family = "Nunito Sans", size = 18,
                                        margin = margin(t = 12, b = 12)),
            legend.position = "right") +
      scale_fill_ddh_d(palette = "protein")

    return(plot_complete)

  }

  #error handling
  tryCatch(make_protein_domain_plot_raw(),
           error = function(x){make_bomb_plot()})
}

# #figure legend
plot_protein_domains_title <- "Protein Domain Plot."
plot_protein_domains_legend <- "Rectangles represent the locations and size of named protein domains, while black shaped elements represent PTMs. Horizontal line(s) indicate the length of one or more selected proteins."

# testing
# make_protein_domain_plot(input = list(content = "ROCK2"), dom_var = "Protein kinase", ptm_var = "N-acetylserine")
# make_protein_domain_plot(input = list(content = c("ROCK1", "ROCK2")), dom_var = "Protein kinase", ptm_var = "N-acetylserine")

## RADIAL PLOT -------------------------------------------------------------
#' Radial Plot
#'
#' \code{make_radial} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a radial plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a radial plot. If an error is thrown, then will return a bomb plot.
#'
#' @export
#' @examples
#' make_radial(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_radial(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_radial(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_radial <- function(cluster_data = sequence_clusters,
                        signature_data = signatures,
                        input = list(),
                        relative = TRUE,
                        cluster = FALSE,
                        barplot = FALSE,
                        card = FALSE) {
  make_radial_raw <- function() {
    # sequence_data_clean <- cluster_data %>%
    #   dplyr::filter(clust != 0)

    sequence_data_clean <- cluster_data

    if(cluster) {
      query_clust <-
        sequence_data_clean %>%
        dplyr::filter(gene_name %in% input$content) %>%
        dplyr::pull(clust) %>%
        unique()

      signature_cluster_means_prep <-
        signature_data %>%
        dplyr::select(uniprot_id, A:Y) %>%
        dplyr::inner_join(cluster_data, by = "uniprot_id") %>%
        dplyr::select(-gene_name, -X1, -X2, -member_prob)

      signature_cluster_means_query <-
        signature_cluster_means_prep %>%
        mutate(clust = as.numeric(as.character(clust)),
               clust = as.factor(ifelse(clust %in% as.numeric(as.character(query_clust)), clust, "Mean"))) %>%
        group_by(clust) %>%
        summarise_if(is.numeric, list(mean = mean)) %>%
        pivot_longer(cols = -clust) %>%
        mutate(name = str_remove(name, "_mean"),
               clust = as.factor(ifelse(clust == "Mean", "Mean", paste0("Cluster ", as.numeric(as.character(clust)))))
        )

      if(relative){
        signature_cluster_means_query <-
          signature_cluster_means_query %>%
          pivot_wider(id_cols = clust, values_fn = mean) %>%
          dplyr::arrange(factor(clust, levels = "Mean")) %>%
          dplyr::mutate(across(A:Y, ~ . / .[1])) %>%
          pivot_longer(cols = -clust)
      }

      if(length(unique(signature_cluster_means_query$clust)) == 1) {
        stop("Unable to cluster this protein by its amino acid sequence so far...")
      }

    } else {
      signature_cluster_means_prep <-
        signature_data %>%
        dplyr::select(gene_name, A:Y)

      signature_cluster_means_query <-
        signature_cluster_means_prep %>%
        dplyr::filter(!gene_name %in% input$content) %>%
        dplyr::rename(Mean = gene_name) %>%
        summarise_if(is.numeric, list(mean = mean)) %>%
        tibble::rownames_to_column("clust") %>%
        pivot_longer(cols = -clust) %>%
        mutate(name = str_remove(name, "_mean"),
               clust = "Mean")

      signature_gene_query <- signature_cluster_means_prep %>%
        dplyr::filter(gene_name %in% input$content) %>%
        pivot_longer(cols = -gene_name) %>%
        dplyr::rename(clust = gene_name)

      signature_cluster_means_query <- dplyr::bind_rows(signature_cluster_means_query,
                                                        signature_gene_query)

      if(relative){
        signature_cluster_means_query <-
          signature_cluster_means_query %>%
          pivot_wider(id_cols = clust, values_fn = mean) %>%
          dplyr::arrange(factor(clust, levels = "Mean")) %>%
          dplyr::mutate(across(A:Y, ~ . / .[1])) %>%
          pivot_longer(cols = -clust)
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
    plot_complete <- ggplot(signature_cluster_means_query,
                            aes(x = fct_inorder(name),
                                y = value,
                                group = clust,
                                color = clust
                            )) +
      {if(!barplot)geom_point(alpha = 0.8, show.legend = FALSE)} +
      {if(!barplot)geom_polygon(fill = NA)} +
      {if(barplot & relative)geom_col(data = signature_cluster_means_query %>%
                                        dplyr::filter(clust != "Mean"),
                                      aes(x = reorder(name, -value),
                                          y = value,
                                          group = clust,
                                          color = clust,
                                          fill = clust),
                                      position = "dodge2")} +
      {if(barplot & !relative)geom_col(aes(fill = clust), position = "dodge2")} +
      {if(barplot & relative)geom_hline(yintercept = 1, color = "gray48")} +
      labs(y = y_label,
           x = element_blank()) +
      {if(!barplot)coord_polar()} +
      scale_color_manual(values = colors_radial) +
      {if(barplot)scale_fill_manual(values = colors_bar)} +
      ## theme changes
      theme_ddh(base_size = 16) +
      {if(!barplot)theme_minimal()} +
      {if(!barplot)theme(
        text = element_text(family = "Nunito Sans"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(family = "Roboto Slab", size = 16),
        axis.text.y = element_text(size = 16, color = "grey30"),
        axis.title = element_text(size = 16)
      )} +
      {if(barplot)theme(
        text = element_text(family = "Nunito Sans"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.text = element_text(family = "Roboto Slab"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
      )} +
      NULL

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        labs(x = "", y = "") + #, title = "Signature Information", caption = "more ...") +
        theme(
          text = element_text(family = "Nunito Sans"),
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_blank(),
          axis.text.x = element_text(family = "Roboto Slab"),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank()
        )
    }

    return(plot_complete)
  }

  #error handling
  tryCatch(make_radial_raw(),
           error = function(x){make_bomb_plot()})
}

#test
# make_radial(input = list(content = "ROCK2"), cluster = F, relative = F, barplot = F, card = F) # default
# make_radial(input = list(content = "ROCK2"), card = TRUE)

#figure legend
plot_radial_title <- "Radial Plot."
plot_radial_legend <- "Amino acid signature/s (percentage of each amino acid in a protein) of the queried gene/s versus the mean amino acid signature of all the other proteins in the dataset (N = 20375)."

plot_aa_bar_title <- "Bar Plot."
plot_aa_bar_legend <- "Amino acid signature/s (percentage of each amino acid in a protein) of the queried gene/s versus the mean amino acid signature of all the other proteins in the dataset (N = 20375)."

cluster_plot_radial_title <- "Cluster Radial Plot."
cluster_plot_radial_legend <- "Amino acid signature/s (percentage of each amino acid in a protein) of the queried cluster/s versus the mean amino acid signature of all the other proteins in the dataset (N = 20375)."

cluster_plot_aa_bar_title <- "Cluster Bar Plot."
cluster_plot_aa_bar_legend <- "Amino acid signature/s (percentage of each amino acid in a protein) of the queried cluster/s versus the mean amino acid signature of all the other proteins in the dataset (N = 20375)."

## UMAP PLOT --------------------------------------------------------
#' UMAP Plot
#'
#' \code{make_umap_plot} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a umap plot plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a umap plot plot. If an error is thrown, then will return a bomb plot.
#'
#' @export
#' @examples
#' make_umap_plot(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_umap_plot(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_umap_plot(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_umap_plot <- function(cluster_data = sequence_clusters,
                           input = list(),
                           show_subset = FALSE,
                           labels = FALSE) {
  make_umap_plot_raw <- function() {

    sequence_data_clean <- cluster_data %>%
      # dplyr::filter(clust != 0) %>%
      mutate(clust = paste0("Cluster ", clust))

    query_clust <- sequence_data_clean %>%
      dplyr::filter(gene_name %in% input$content) %>%
      dplyr::pull(clust) %>%
      unique()

    cluster_genes <- sequence_data_clean %>%
      dplyr::filter(clust %in% query_clust)

    colors <- ddh_pal_c(palette = "protein")(length(query_clust))

    # UMAP PLOT
    plot_complete <- ggplot() +
      {if(!show_subset)geom_point(data = sequence_data_clean %>%
                                    dplyr::filter(!uniprot_id %in% cluster_genes$uniprot_id),
                                  aes(X1, X2), size = 0.8, color = "grey80")} +
      geom_point(data = sequence_data_clean %>%
                   dplyr::filter(uniprot_id %in% cluster_genes$uniprot_id),
                 aes(X1, X2, color = clust), size = 0.8) +
      {if(labels)ggrepel::geom_label_repel(data = cluster_genes %>%
                                             dplyr::filter(gene_name %in% input$content),
                                           aes(X1, X2, label = gene_name))} +
      labs(x = "UMAP1",
           y = "UMAP2") +
      scale_color_manual(
        values = rep(colors, length.out =
                       nrow(sequence_data_clean %>%
                              dplyr::filter(uniprot_id %in% cluster_genes$uniprot_id)))
      ) +
      guides(color = guide_legend(override.aes = list(size = 3))) +
      ## theme changes
      theme_ddh() +
      theme(
        text = element_text(family = "Nunito Sans"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.text = element_text(family = "Roboto Slab", size = 16),
        axis.text.y = element_text(size = 16, color = "grey30"),
        axis.title = element_text(size = 16)
      ) +
      NULL

    return(plot_complete)
  }

  #error handling
  tryCatch(make_umap_plot_raw(),
           error = function(x){make_bomb_plot()})
}

#figure legend
plot_umap_title <- "UMAP Plot."
plot_umap_legend <- "Amino acid signature (percentage of each amino acid in a protein) UMAP embeddings (2D) colored by the cluster to which they belong (N = 20375)."

## CLUSTER ENRICHMENT PLOT --------------------------------------------------------
#' Cluster Enrich Plot
#'
#' \code{make_cluster_enrich} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a cluster enrich plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a cluster enrich plot. If an error is thrown, then will return a bomb plot.
#'
#' @export
#' @examples
#' make_cluster_enrich(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_cluster_enrich(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_cluster_enrich(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_cluster_enrich <- function(input = list(),
                                ontology = "BP",
                                num_terms = 20) {

  make_cluster_enrich_plot_raw <- function() {

    plot_data <- make_clustering_enrichment_table(input = input, #update this
                                                  ontology = ontology)

    plot_complete <-
      plot_data %>%
      dplyr::arrange(pvalue) %>%
      dplyr::slice(1:num_terms) %>%
      ggplot(aes(x = Count,
                 y = reorder(substr(paste0(Description, " (", ID, ")"), 1, 40), Count),
                 fill = pvalue)) +
      geom_col() +
      labs(x = "Gene Count",
           y = NULL) +
      scale_fill_ddh_c(palette = "protein", reverse = TRUE) +
      guides(fill = guide_colorbar(barheight = unit(6, "lines"),
                                   barwidth = unit(.6, "lines"),
                                   reverse = TRUE)) +
      theme_ddh() +
      theme(
        title = element_blank(),
        text = element_text(family = "Nunito Sans"),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.text = element_text(family = "Roboto Slab"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
      ) +
      NULL

    return(plot_complete)
  }
  #error handling
  tryCatch(make_cluster_enrich_plot_raw(),
           error = function(x){make_bomb_plot()})
}

# #figure legend
plot_cluster_enrichment_title <- "Cluster Enrichment Plot."
plot_cluster_enrichment_legend <- "Enriched GO terms (BP, MF, and CC) by all genes in the selected amino acid signature cluster. The x-axis shows the number of genes in the cluster that belong to each term while the color scale represents the p-values of each enriched term."

## STRUCTURE PLOT --------------------------------------------------------------------
#' Structure Plot
#'
#' \code{make_structure} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a structure plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a structure plot. If an error is thrown, then will return a bomb plot.
#'
#' @export
#' @examples
#' make_structure(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_structure(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_structure(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_structure <- function(input = list(),
                           card = FALSE) {
  make_structure_raw <- function() {
    #if multiple, then pull single "rep" image; consider pulling >1 and using patchwork, eg.
    if(length(input$content > 1)){
      fav_gene <- input$content[[1]] #sample(input$content, 1)
    } else {
      fav_gene <- input$content
    }

    #fetch image from dir
    file_name <- glue::glue('{fav_gene}_structure_plot.jpg')
    image_path <- here::here("data", "images", "gene", fav_gene, file_name)

    if(card == TRUE&&file.exists(image_path)){
      protein_image <-
        magick::image_read(path = image_path)

      #scale up if too small to make a square
      image_details <- magick::image_info(protein_image)

      #center when scaling by providing offset
      orig_width <- image_details$width
      orig_height <- image_details$height
      #target
      width <- 1080
      height <- 1080
      #split difference in pixels between scale
      offset_x <- as.integer((width - orig_width) / 2)
      offset_y <- as.integer((height - orig_height) / 2)
      #build offsets
      width_offset <- glue::glue('{width}x+{offset_x}')
      height_offset <- glue::glue('x{height}+{offset_y}')
      #perform scaling
      if(image_details$width < 1080) {protein_image <- image_scale(protein_image, width_offset)}
      if(image_details$height < 1080) {protein_image <- image_scale(protein_image, height_offset)}

      protein_image <-
        protein_image %>%
        image_crop("1080x1080")

      #place image on ggplot
      plot_complete <-
        ggplot() +
        ggpubr::background_image(protein_image) +
        labs(x = "") + #, title = "Structure Information", caption = "more ...") +
        theme_void()

      return(plot_complete)
    }

    if(file.exists(image_path)){return(image_path)} else {return(NULL)}

  }
  #error handling
  tryCatch(make_structure_raw(),
           error = function(x){make_bomb_plot()})
}

#make_structure(input = list(content = "ROCK1"))
#make_structure(input = list(content = "ROCK1"), card = TRUE)

protein_structure_title <- "Predicted Structure."
protein_structure_legend <- "Alpha Fold predicted structure rendered in a ribbon diagram."

## 3D STRUCTURE PLOT --------------------------------------------------------------------
make_structure3d <- function(pdb_ids = uniprot_pdb_table,
                             protein_data = proteins,
                             gene_id = NULL,
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

    if(is.null(gene_id)) {
      if(length(input$content) == 1) {
        gene_id <- input$content
      } else {
        gene_id <- input$content[1]
      }
    }

    pdb_path <- load_pdb(input = list(content = gene_id))

    if(!is.null(pdb_path) & is.null(pdb_id)) {
      plot_data <- read.pdb(pdb_path)
    } else {
      plot_data <- protein_data %>%
        filter(gene_name %in% gene_id) %>%
        left_join(pdb_ids, by = c("uniprot_id" = "uniprot")) %>%
        unnest(data) %>%
        {if (is.null(pdb_id)) dplyr::slice(., 1) else filter(., pdb == pdb_id)} %>%
        pull(pdb) %>%
        r3dmol::m_fetch_pdb()
    }

    plot_complete <- r3dmol::r3dmol() %>%
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
        m_add_style(
          style = c(
            m_style_stick(),
            m_style_sphere(scale = 0.3)
          ),
          sel = m_sel(resi =  eval(parse(text = resi)),
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
           error = function(x){make_bomb_plot()})
}

# make_structure3d(input = list(content = "ROCK1"))

protein_structure_title3d <- "Predicted 3D Structure."
protein_structure_legend3d <- "Interactive 3D predicted structure."

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
#' @export
#' @examples
#' make_pubmed(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_pubmed(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_pubmed(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_pubmed <- function(pubmed_data = pubmed,
                        input = list(),
                        card = FALSE) {
  make_pubmed_raw <- function() {
    plot_data <-
      pubmed_data %>%
      dplyr::filter(name %in% input$content) %>%
      tidyr::unnest(data) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(name, year) %>%
      dplyr::summarize(n = n()) %>%
      dplyr::mutate(cumsum = cumsum(n)) %>%
      dplyr::arrange(year)

    plot_max <-
      plot_data %>%
      dplyr::filter(!is.na(year)) %>%
      dplyr::filter(cumsum == max(cumsum))

    plot_step <-
      ggplot(plot_data) +
      geom_step(
        aes(x = year,
            y = cumsum,
            group = name,
            color = fct_reorder2(name, year, cumsum)),
        size = 1.2
      ) +
      coord_cartesian(clip = "off") + #allows points & labels to fall off plotting area
      scale_x_continuous(breaks = scales::breaks_pretty(4),
                         expand = c(0, 0)) +
      scale_y_continuous(expand = c(.01, .01), limits = c(0, max(plot_max$cumsum))) +
      scale_color_ddh_d(palette = input$type) +
      labs(x = "Year of Publication", y = "Cumulative Sum") +
      theme_ddh(base_size = 16,
                margin = 20) +
      theme(
        #text = element_text(family = "Nunito Sans"),
        #axis.text.y = element_text(family = "Roboto Slab"),
        panel.grid.major.y = element_line(color = "grey75", size = .5, linetype = "15")
      ) +
      NULL

    ## only add labels + white space if several genes queried
    if (length(unique(plot_data$name)) > 1) {
      plot_complete <-
        plot_step +
        {if(card == FALSE)
          geom_text_repel(
            data = plot_max,
            aes(x = year,
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
        geom_point(
          data = plot_max,
          aes(x = year,
              y = cumsum,
              group = name,
              color = fct_reorder2(name, year, cumsum)),
          size = 4, shape = 21, fill = "white", stroke = 2
        ) +
        theme(
          plot.margin = margin(7, 180, 7, 7) #adds margin to right side of graph for label #adds margin to right side of graph for label
        )
    } else {
      plot_complete <-
        plot_step +
        geom_point(
          data = plot_max,
          aes(x = year,
              y = cumsum,
              group = name,
              color = fct_reorder2(name, year, cumsum)),
          size = 4, shape = 21, fill = "white", stroke = 2
        )
    }

    if(length(input$content) == 1){ #Fix this when data()$content is updated (need to be generic to handle multiple data types)
      plot_complete  <-
        plot_complete +
        labs(y = "Cumulative Publications",
             color = "") +
        guides(color = "none")
    } else {
      plot_complete <-
        plot_complete +
        labs(y = "Cumulative Publications",
             color = "Query")
    }

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        theme(plot.margin = margin(5, 10, 5, 5),
              legend.position="none") +
        labs(x = "") +
        NULL
    }

    return(plot_complete)
  }
  #error handling
  tryCatch(make_pubmed_raw(),
           error = function(x){make_bomb_plot()})
}

# make_pubmed(input = list(type = "gene", query = "ROCK1", content = "ROCK1"))
# make_pubmed(input = list(type = "gene", query = "ROCK1", content = "ROCK1"), card = TRUE)
# make_pubmed(input = list(type = "gene", query = "custom_gene_list", content = c("ROCK1", "ROCK2")), card = TRUE)
# make_pubmed(input = list(type = "gene", query = "custom_gene_list", content = c("RDX", "ROCK2", "DTX3L", "MSN", "SORL1", "EZR")), card = TRUE)
# make_pubmed(input = list(type = "pathway", query = "0060148", content = c("FXR1", "ZFP36", "DHX9", "XPO5", "FMR1", "STAT3", "WTIP", "PUM2", "AJUBA", "PUM1", "LIMD1")))
# make_pubmed(input = list(type = "compound", content = "aspirin"))
# make_pubmed(input = list(type = "gene", query = "ROCK1", content = "ROCK1"), card = TRUE)

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
#' @export
#' @examples
#' make_cellanatogram(input = list(type = "gene", content = c("ROCK2")))
#' make_cellanatogram(input = list(type = "gene", content = c("ROCK2")), card = TRUE)
#' make_cellanatogram(input = list(type = "gene", query = "ROCK1", content = c("ROCK1", "ROCK2")))
#' \dontrun{
#' make_cellanatogram(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_cellanatogram <- function(cellanatogram_data = subcell,
                               input = list(),
                               card = FALSE) {
  make_cellanatogram_raw <- function() {
    plot_data <-
      cellanatogram_data %>%
      dplyr::filter_all(any_vars(gene_name %in% input$content)) %>%
      dplyr::filter(!is.na(type)) %>%
      dplyr::select(-value) %>%
      dplyr::add_count(main_location) %>%
      dplyr::mutate(value = as_factor(n)) %>%
      #dplyr::arrange(desc(row_number())) ## not sure this always works...
      dplyr::mutate(
        organ = stringr::str_replace_na(organ, "cytosol"),
        type = factor(gene_name),
        organ = fct_rev(organ)
      )

    plot_complete <-
      plot_data %>%
      gganatogram(outline = TRUE, fillOutline = 'grey95', organism = "cell", fill = "value") +
      theme_void(base_size = 14) +
      theme(
        text = element_text(family = "Nunito Sans"),
        plot.margin = margin(5, 10, 5, 5)
      ) +
      coord_fixed() +
      scale_fill_ddh_d(palette = "protein") +
      labs(fill = "Count") +
      NULL

    if(length(input$content) == 1){
      plot_complete  <-
        plot_complete +
        guides(fill = "none")
    }

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        labs(x = "") +
        guides(fill = "none") +
        NULL
    }
    return(plot_complete)
  }
  #error handling
  tryCatch(make_cellanatogram_raw(),
           error = function(x){make_bomb_plot()})
}

# make_cellanatogram(input = list(type = "gene", content = c("ROCK2")))
# make_cellanatogram(input = list(type = "gene", content = c("ROCK2")), card = TRUE)
# make_cellanatogram(input = list(type = "gene", query = "ROCK1", content = c("ROCK1", "ROCK2")))


# make cell anatogram facet
make_cellanatogramfacet <- function(cellanatogram_data = subcell,
                                    input = list()) {
  make_cellanatogramfacet_raw <- function() {
    plot_complete <-
      cellanatogram_data %>%
      dplyr::filter(gene_name %in% input$content,
                    !is.na(type)) %>%
      dplyr::mutate(colour = scales::alpha(ddh_pal_d(palette = "protein")(1), .75)) %>%
      full_join(tibble(gene_name = input$content)) %>%
      dplyr::mutate(
        #colour = if_else(is.na(organ), "grey92", scales::alpha(ddh_pal_d(palette = "gene")(1), .75)),
        organ = stringr::str_replace_na(organ, "cytosol"),
        type = factor(gene_name),
        organ = fct_rev(organ)
      ) %>%
      gganatogram(outline = TRUE, fillOutline = 'grey95', organism = "cell", fill = "colour") +
      theme_void(base_size = 14) +
      theme(
        text = element_text(family = "Nunito Sans", face = "bold"),
        plot.margin = margin(5, 10, 5, 5),
        panel.spacing.y = unit(1.2, "lines")
      ) +
      facet_wrap(~ type, ncol = 3, drop = FALSE) +
      coord_fixed() +
      labs(fill = "Count") +
      theme(panel.spacing.y = unit(1.2, "lines"),
            strip.text = element_text(size = 18, family = "Roboto Slab", face = "plain")) +
      NULL

    return(plot_complete)
  }
  #error handling
  tryCatch(make_cellanatogramfacet_raw(),
           error = function(x){make_bomb_plot()})
}

# make_cellanatogramfacet(input = list(type = "gene", content = c("ROCK1", "ROCK2")))

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
#' @export
#' @examples
#' make_female_anatogram(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_female_anatogram(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_female_anatogram(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_female_anatogram <- function(anatogram = "female",
                                  input = list(),
                                  card = FALSE) {
  make_female_anatogram_raw <- function() {
    if(anatogram == "female"){
      raw_body_data = female_tissue
    } else if(anatogram == "male") {
      raw_body_data = male_tissue
    } else {
      return("declare your anatogram type")
    }
    body_data <-
      raw_body_data %>%
      dplyr::filter_all(any_vars(gene_name %in% input$content)) %>%
      dplyr::filter(!is.na(type)) %>%
      dplyr::arrange(dplyr::desc(-value))

    break_points <- unname(quantile(body_data$value))

    if(length(input$content) > 1){
      body_data <-
        body_data %>%
        group_by(organ) %>%
        mutate(value=sum(value)) %>%
        arrange(desc(-value))
    }

    body_plot <-
      body_data %>%
      gganatogram(outline = TRUE, fillOutline='grey95', organism = "human", sex = anatogram, fill = 'value') +
      coord_fixed() +
      scale_fill_ddh_c(palette = "gene", breaks = break_points, name = NULL) +
      guides(fill = guide_colorsteps(show.limits = TRUE)) +
      theme_void(base_size = 16, base_family = "Chivo") +
      NULL

    if(card == TRUE){
      body_plot <-
        body_plot +
        labs(x = "") + #, title = "Tissue Distribution", caption = "more ...") +
        NULL
    }

    return(body_plot)
  }
  #error handling
  tryCatch(make_female_anatogram_raw(),
           error = function(x){make_bomb_plot()})
}

# make_female_anatogram(input = list(type = "gene", content = c("ROCK1")))
# make_female_anatogram(input = list(type = "gene", content = c("ROCK1")), card = TRUE)
# make_female_anatogram(input = list(type = "gene", query = "ROCK1", content = c("ROCK1", "ROCK2")))


make_male_anatogram <- function(anatogram = "male",
                                input = list(),
                                card = FALSE){
  male_anatogram <- make_female_anatogram(anatogram,
                                          input,
                                          card)
  return(male_anatogram)
}

# make_male_anatogram(input = list(type = "gene", content = c("ROCK1")))
# make_male_anatogram(input = list(type = "gene", content = c("ROCK1")), card = TRUE)

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
#' @export
#' @examples
#' make_tissue(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_tissue(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_tissue(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_tissue <- function(tissue_data = tissue,
                        input = list(),
                        card = FALSE) {
  make_tissue_raw <- function() {
    plot_data <-
      tissue_data %>%
      dplyr::filter_all(any_vars(gene_name %in% input$content)) %>%
      group_by(organ) %>%
      dplyr::mutate(sum_value = sum(value),
                    organ = stringr::str_replace_all(organ, "_", " "),
                    organ = stringr::str_to_title(organ)) %>%
      arrange(desc(sum_value))

    if(nrow(plot_data) == 0){return(NULL)}

    if(card == TRUE){
      #get distinct
      top_organs <- plot_data %>%
        distinct(organ) %>%
        ungroup() %>%
        slice_head(n = 10) %>%
        pull()
      plot_data <-
        plot_data %>%
        dplyr::filter(organ %in% top_organs)
    }

    plot_draft <-
      ggplot(plot_data,
             aes(x = value,
                 y = fct_reorder(organ, sum_value))) +
      coord_cartesian(clip = "off") +
      scale_x_continuous(expand = c(0, 0), sec.axis = dup_axis()) +
      scale_y_discrete(expand = c(.01, .01)) +
      theme_ddh() +
      theme(axis.text.y = element_text(angle = 0, family = "Roboto Slab"),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank()) +
      labs(y = NULL)

    ## Version with color mapped to count for single gene queries
    if(length(input$content) == 1){
      plot_complete <-
        plot_draft +
        geom_col(aes(fill = value), width = .82) +
        scale_fill_ddh_c(palette = "gene", guide = "none") +
        labs(x = paste0(str_c(input$content, collapse = ", "), " Normalized Expression")) +
        theme(plot.margin = margin(0, 15, 0, 0)) # add some space to avoid cutting of labels
    } else {
      plot_complete <-
        plot_draft +
        geom_col(aes(fill = gene_name), width = .82) +
        scale_fill_ddh_d(palette = "gene", shuffle = TRUE, seed = 5L) +
        labs(x = "Sum of Normalized Expression",
             fill = "Query\nGene") +
        theme(legend.justification = "top")
    }

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        labs(x = "") +
        theme(axis.text.x=element_blank(),
              legend.position='none') +
        NULL
    }

    return(plot_complete)
  }
  #error handling
  tryCatch(make_tissue_raw(),
           error = function(x){make_bomb_plot()})
}

# make_tissue(input = list(type = "gene", content = c("ROCK1")))
# make_tissue(input = list(type = "gene", content = c("ROCK1")), card = TRUE)
# make_tissue(input = list(type = "gene", content = c("ROCK1", "ROCK2")))
# make_tissue(input = list(type = "gene", content = c("ROCK1", "ROCK2")), card = TRUE)
# make_tissue(input = list(type = "gene", content = c("RDX", "ROCK2", "DTX3L", "MSN", "SORL1", "EZR")), card = TRUE)




## EXPRESSION PLOT --------------------------------------------------------
#' Cell Expression Plot
#'
#' \code{make_cellexpression} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a cell expression "rug" plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a cellexpression plot. If an error is thrown, then will return a bomb plot.
#'
#' @export
#' @examples
#' make_cellexpression(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_cellexpression(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_cellexpression(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_cellexpression <- function(expression_data = expression_long,
                                expression_join = expression_names,
                                input = list(),
                                var = "gene",
                                card = FALSE) {
  make_cellexpression_raw <- function() {
    if (var == "gene") {
      plot_initial <-
        expression_data %>%
        dplyr::select(any_of(c("X1", "gene", "gene_expression"))) %>%
        dplyr::rename("expression_var" = "gene_expression")
      mean <- mean_virtual_gene_expression
      upper_limit <- gene_expression_upper
      lower_limit <- gene_expression_lower
    } else if (var == "protein") {
      plot_initial <-
        expression_data %>%
        dplyr::select(X1, gene, protein_expression) %>%
        dplyr::rename("expression_var" = "protein_expression")
      mean <- mean_virtual_protein_expression
      upper_limit <- protein_expression_upper
      lower_limit <- protein_expression_lower
    } else {
      stop("declare your variable")
    }

    if (input$type == "gene") {
      plot_data <- plot_initial %>%
        dplyr::filter(gene %in% input$content,
                      !is.na(expression_var)) %>%
        dplyr::left_join(expression_join, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::select(cell_line, lineage, lineage_subtype, everything()) %>%
        dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
        #dplyr::mutate(gene_fct = fct_reorder(gene, expression_var, .fun = max, .desc = TRUE)) %>%
        dplyr::mutate(gene_fct = fct_inorder(gene)) %>%
        ggplot(aes(y = gene_fct,
                   x = expression_var,
                   text = paste0("Cell Line: ", cell_line),
                   color = gene
        ))
    } else if (input$type == "cell") {
      plot_data <- plot_initial %>%
        dplyr::left_join(expression_join, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::select(cell_line, lineage, lineage_subtype, everything()) %>%
        dplyr::filter(cell_line %in% input$content,
                      !is.na(expression_var)) %>%
        dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
        dplyr::mutate(cell_fct = fct_inorder(cell_line)) %>%
        ggplot(aes(y = cell_fct,
                   x = expression_var,
                   text = paste0("Gene: ", gene),
                   color = cell_line
        ))
    }

    plot_complete <-
      plot_data +
      geom_point(alpha = 0.1, shape = "|", size = 12) +
      geom_vline(xintercept = 0, color = "lightgray") + #3SD is below zero
      geom_vline(xintercept = mean) +
      geom_vline(xintercept = upper_limit, color = "lightgray") +#3SD
      scale_x_continuous(expand = expansion(mult = 0.01)) +
      scale_y_discrete(expand = expansion(mult = 1 / length(input$content)), na.translate = FALSE) +
      scale_color_ddh_d(palette = input$type) +
      theme_ddh() +
      theme(
        text = element_text(family = "Nunito Sans"),
        axis.text = element_text(family = "Roboto Slab"),
        axis.text.y = element_text(size = 18),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
      ) +
      NULL

    if(length(input$content) == 1){
      plot_complete  <-
        plot_complete +
        theme(axis.text.y = element_blank()) +
        labs(x = paste0(str_c(input$content, collapse = ", "), " ", var, " levels"))
    } else {
      plot_complete <-
        plot_complete +
        labs(x = paste0(var, " levels"), color = "Query")
    }

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        labs(x = "")
    }
    return(plot_complete)
  }
  #error handling
  tryCatch(print(make_cellexpression_raw()),
           error = function(x){make_bomb_plot()})
}

# make_cellexpression(input = list(type = "gene", content = c("ROCK1")))
# make_cellexpression(input = list(type = "gene", content = c("ROCK1")), var = "protein")
# make_cellexpression(input = list(type = "gene", content = c("ROCK1")), var = "protein")
# make_cellexpression(input = list(type = "gene", content = c("ROCK1", "ROCK2")))
# make_tissue(input = list(type = "gene", content = c("ROCK1")), card = TRUE)
# make_cellexpression(input = list(type = "gene", content = c("HSP90B1")))


#figure legend
plot_cellexp_title <- "Expression Values."
#plot_cellexp_legend <- paste0("Each point shows the ranked expression value across ", n_distinct(expression_long$X1)," cell lines.")
plot_cellexp_legend <- "Each point shows the ranked expression value across CCLE cell lines."
plot_cellLineexp_title <- "Expression Values."
#plot_cellLineexp_legend <- paste0("Each point shows the ranked expression value across ", n_distinct(expression_long$gene)," genes.")
plot_cellLineexp_legend <- "Each point shows the ranked expression value across genes."
#If present, black line indicates resampled mean expression value (", round(mean_virtual_expression, 2), "). If present, gray line indicates 3 standard deviations away from the resampled mean (", round(expression_upper, 2), ").

# G-EXPvP-EXP ----------------------------------
#this makes gene v. protein in a cell
make_cellgeneprotein <- function(expression_data = expression_long,
                                 expression_join = expression_names,
                                 input = list(),
                                 card = FALSE) {
  make_cellgeneprotein_raw <- function() {
    if (input$type == "gene") {
      plot_initial <- expression_data %>%
        dplyr::filter(gene %in% input$content,
                      !is.na(gene_expression),
                      !is.na(protein_expression)) %>%
        dplyr::left_join(expression_join, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::select(cell_line, lineage, lineage_subtype, everything()) %>%
        dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
        ggplot(aes(x = gene_expression,
                   y = protein_expression,
                   text = paste0("Cell Line: ", cell_line),
                   color = gene,
                   group = gene)
        )
    } else if (input$type == "cell") {
      plot_initial <- expression_data %>%
        dplyr::left_join(expression_join, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::select(cell_line, lineage, lineage_subtype, everything()) %>%
        dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
        dplyr::filter(cell_line %in% input$content,
                      !is.na(gene_expression),
                      !is.na(protein_expression)) %>%
        ggplot(aes(x = gene_expression,
                   y = protein_expression,
                   text = paste0("Gene: ", gene),
                   color = cell_line,
                   group = cell_line)
        )
    }

    plot_complete <-
      plot_initial +
      geom_point(alpha = 0.4) +
      #add geom to drop linear regression line?
      geom_hline(yintercept = 0, color = "lightgray") +
      geom_vline(xintercept = 0, color = "lightgray") +
      # smooth line
      geom_smooth(method = "lm",
                  se = TRUE) +
      # R coefs
      {if(card == FALSE)ggpubr::stat_cor(digits = 3)} +
      scale_color_ddh_d(palette = input$type) +
      theme_ddh() +
      theme(
        text = element_text(family = "Nunito Sans"),
        axis.text = element_text(family = "Roboto Slab")
      ) +
      NULL

    if(length(input$content) == 1){
      plot_complete  <-
        plot_complete +
        labs(x = paste0(input$content, " Gene Expression"), y = paste0(input$content, " Protein Expression")) +
        theme(legend.position = "none")
    } else if(input$type == "gene") {
      plot_complete <-
        plot_complete +
        labs(x = "Gene Expression", y = "Protein Expression", color = "Query \nGene")
    } else if(input$type == "cell") {
      plot_complete <-
        plot_complete +
        labs(x = "Gene Expression", y = "Protein Expression", color = "Query \nCell Line")
    }

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        labs(x = "", y = "") +
        theme(legend.position='none')
    }
    return(plot_complete)
  }
  #error handling
  tryCatch(print(make_cellgeneprotein_raw()),
           error = function(x){make_bomb_plot()})
}

plot_cellgeneprotein_title <- "Gene Expression versus Protein Expression."
plot_cellgeneprotein_legend <- "Each point shows the gene expression value compared to the protein expression value for gene within a given cell line. The Pearson correlation coefficient and the p-values are provided in the top-left corner of the plot."

## CELL DEPS --------------------------------------------------------------------
#' Cell Dependencies Plot
#'
#' \code{make_celldeps} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a celldeps plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a celldeps plot. If an error is thrown, then will return a bomb plot.
#'
#' @export
#' @examples
#' make_celldeps(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_celldeps(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_celldeps(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_celldeps <- function(celldeps_data = achilles_long,
                          prism_data = prism_long,
                          expression_data = expression_names,
                          input = list(),
                          card = FALSE,
                          lineplot = FALSE,
                          scale = NULL) {#scale is expecting 0 to 1
  make_celldeps_raw <- function() {
    if(input$type == "gene") {
      aes_var <- rlang::sym("name")
      var_title <- "Gene"
      ylab <- "Dependecy Score"
      mean <- mean_virtual_achilles

      plot_data <-
        celldeps_data %>% #plot setup
        dplyr::filter(gene %in% input$content) %>%
        dplyr::left_join(expression_data, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::group_by(gene) %>%
        dplyr::arrange(dep_score) %>%
        dplyr:: mutate(
          rank = 1:n(),
          med = median(dep_score, na.rm= TRUE)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::rename(name = gene)

    } else if(input$type == "compound") {
      aes_var <- rlang::sym("name")
      var_title <- "Compound"
      ylab <- "Log2FC"
      mean <- mean_virtual_achilles_cell_line
      plot_data <-
        prism_data %>% #plot setup
        dplyr::filter(name %in% input$content) %>%
        dplyr::left_join(expression_data, by = c("x1" = "X1")) %>%
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
      mean <- mean_virtual_achilles_cell_line

      plot_data <-
        celldeps_data %>% #plot setup
        dplyr::left_join(expression_data, by = "X1") %>%
        dplyr::filter(cell_line %in% input$content) %>%
        dplyr::select(-X1) %>%
        dplyr::group_by(cell_line) %>%
        dplyr::arrange(dep_score) %>%
        dplyr:: mutate(
          rank = 1:n(),
          med = median(dep_score, na.rm= TRUE)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::rename(name = gene)
    }

    if(!is.null(scale) & !card){
      plot_data <-
        plot_data %>%
        slice_sample(prop = scale)
    }

    if(card) {
      if(is.null(scale)) {
        scale <- 0.3
      }
      if(input$type == "gene") {
        plot_data <- plot_data %>%
          group_by(name) %>%
          sample_n(scale*n()) %>%
          ungroup()
      } else {
        plot_data <- plot_data %>%
          group_by(cell_line) %>%
          sample_n(scale*n()) %>%
          ungroup()
      }
    }

    plot_complete <-
      plot_data %>%
      ggplot(aes(x = rank,
                 y = dep_score,
                 text = glue::glue('{var_title}: {name}\nCell Line: {cell_line}'),
                 color = fct_reorder(!!aes_var, med),
                 fill = fct_reorder(!!aes_var, med)
      )) +
      ## dot/line plot
      {if(!card & !lineplot)geom_point(size = 1.1, stroke = .25, alpha = 0.6)} +
      {if(input$type == "gene" & (card | lineplot))geom_line(aes(group = name))} +
      {if(input$type == "cell" & (card | lineplot))geom_line(aes(group = cell_line))} +
      ## indicator lines dep. score
      geom_hline(yintercept = mean) +
      geom_hline(yintercept = 1, size = .2, color = "grey70", linetype = "dashed") +
      geom_hline(yintercept = -1, size = .2, color = "grey70", linetype = "dashed") +
      geom_hline(yintercept = 0, size = .2, color = "grey50") +
      ## scales + legends
      #scale_x_discrete(expand = expansion(mult = 0.02), na.translate = FALSE) +
      scale_color_ddh_d(palette = input$type) +
      scale_fill_ddh_d(palette = input$type) +
      guides(
        color = guide_legend(reverse = TRUE, override.aes = list(size = 4, stroke = .8)),
        fill = guide_legend(reverse = TRUE, override.aes = list(size = 3.8, stroke = .8))
      ) +
      ## titles
      labs(
        x = NULL,
        y = ylab,
        color = "Query",
        fill = "Query"
      ) +
      ## theme changes
      theme_ddh() +
      theme(
        text = element_text(family = "Nunito Sans"),
        axis.text = element_text(family = "Roboto Slab"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
      ) +
      NULL

    ##only plot legend in case of more than 1 gene selected
    if(length(input$content) == 1){
      plot_complete  <-
        plot_complete +
        theme(legend.position = "none")
      plot_complete
    } else {
      plot_complete
    }

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        labs(x = "") +
        theme(legend.position = "none")
    }

    return(plot_complete)
  }
  #error handling
  tryCatch(make_celldeps_raw(),
           error = function(x){make_bomb_plot()})
}

#figure legend
plot_celldeps_title <- "Dependency Curve."
plot_celldeps_legend <- "Each point shows the ranked dependency score ordered from low to high scores. Dependency scores less than -1 indicate a gene that is essential within a cell line. Dependency scores close to 0 mean no changes in fitness when the gene is knocked out. Dependency scores greater than 1 indicate gene knockouts lead to a gain in fitness."

#test
#make_celldeps(input = list(type = "gene", query = "ROCK1", content = c("ROCK1", "ROCK2")), mean = mean_virtual_achilles)
#make_celldeps(input = list(type = "gene", query = "ROCK1", content = c("ROCK1")), mean = mean_virtual_achilles, card = TRUE)
#make_celldeps(input = list(type = "compound", query = "aspirin", content = "aspirin"), mean = mean_virtual_achilles_cell_line)
#make_celldeps(input = list(type = "cell", query = "HEPG2", content = "HEPG2"), card = TRUE)
#make_celldeps(input = list(type = "cell", query = "HEPG2", content = "HEPG2"), scale = 0.5)
#make_celldeps() #error

# tic()
# make_celldeps(input = list(type = "cell", query = "HEPG2", content = "HEPG2"), scale = 0.01)
# toc()

## BAR PLOT --------------------------------------------------------------------
#' Cell Dependencies Bar Plot
#'
#' \code{make_cellbar} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a cellbar plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a cellbar plot. If an error is thrown, then will return a bomb plot.
#'
#' @export
#' @examples
#' make_cellbar(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_cellbar(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_cellbar(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_cellbar <- function(celldeps_data = achilles_long,
                         prism_data = prism_long,
                         expression_data = expression_names,
                         input = list(),
                         mean,
                         card = FALSE,
                         scale = NULL) {
  make_cellbar_raw <- function() {
    if(input$type == "gene") {
      aes_var <- rlang::sym("name")
      var_title <- "Gene"
      ylab <- "Dependecy Score"
      mean <- mean_virtual_achilles

      plot_data <-
        celldeps_data %>% #plot setup
        dplyr::filter(gene %in% input$content) %>%
        dplyr::left_join(expression_data, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::group_by(gene) %>%
        dplyr::arrange(dep_score) %>%
        dplyr:: mutate(
          rank = 1:n(),
          med = median(dep_score, na.rm= TRUE)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::rename(name = gene) %>%
        mutate(rank = as.integer(fct_reorder(cell_line, dep_score)))

    } else if(input$type == "compound") {
      aes_var <- rlang::sym("name")
      var_title <- "Compound"
      ylab <- "Log2FC"
      mean <- mean_virtual_prism_cor

      plot_data <-
        prism_data %>% #plot setup
        dplyr::filter(name %in% input$content) %>%
        dplyr::left_join(expression_data, by = c("x1" = "X1")) %>%
        dplyr::select(-1) %>%
        dplyr::group_by(name) %>%
        dplyr::arrange(log2fc) %>%
        dplyr:: mutate(
          rank = 1:n(),
          med = median(log2fc, na.rm= TRUE)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::rename(dep_score = log2fc) %>% #rename for graph
        mutate(rank = as.integer(fct_reorder(name, dep_score)))

    } else { #cell lines
      aes_var <- rlang::sym("cell_line")
      var_title <- "Gene"
      ylab <- "Dependecy Score"
      mean <- mean_virtual_achilles_cell_line

      plot_data <-
        celldeps_data %>% #plot setup
        dplyr::left_join(expression_data, by = "X1") %>%
        dplyr::filter(cell_line %in% input$content) %>%
        dplyr::select(-X1) %>%
        dplyr::group_by(cell_line) %>%
        dplyr::arrange(dep_score) %>%
        dplyr:: mutate(
          rank = 1:n(),
          med = median(dep_score, na.rm= TRUE)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::rename(name = gene) %>%
        mutate(rank = as.integer(fct_reorder(name, dep_score)))

    }

    if(card) {
      if(is.null(scale)) {
        scale <- 0.3
      }
      if(input$type == "gene") {
        plot_data <- plot_data %>%
          group_by(name) %>%
          sample_n(scale*n()) %>%
          ungroup()
      } else {
        plot_data <- plot_data %>%
          group_by(cell_line) %>%
          sample_n(scale*n()) %>%
          ungroup()
      }
    }

    plot_complete <-
      plot_data %>%
      ggplot(aes(x = rank,
                 y = dep_score,
                 text = glue::glue('{var_title}: {name}\nCell Line: {cell_line}'),
                 color = fct_reorder(!!aes_var, med),
                 fill = fct_reorder(!!aes_var, med)
      )) +
      ## bar plot
      geom_bar(stat = "identity", width = 0.5) +
      ## indicator lines dep. score
      geom_hline(yintercept = mean) +
      geom_hline(yintercept = 1, size = .2, color = "grey70", linetype = "dashed") +
      geom_hline(yintercept = -1, size = .2, color = "grey70", linetype = "dashed") +
      geom_hline(yintercept = 0, size = .2, color = "grey50") +
      ## scales + legends
      scale_x_discrete(expand = expansion(mult = 0.02), na.translate = FALSE) +
      scale_color_ddh_d(palette = input$type) +
      scale_fill_ddh_d(palette = input$type) +
      guides(
        color = guide_legend(reverse = TRUE, override.aes = list(size = 4, stroke = .8)),
        fill = guide_legend(reverse = TRUE, override.aes = list(size = 3.8, stroke = .8))
      ) +
      ## titles
      labs(
        x = NULL,
        y = ylab,
        color = "Query",
        fill = "Query"
      ) +
      ## theme changes
      theme_ddh() + #base_size = 15 default
      theme(
        text = element_text(family = "Nunito Sans"),
        axis.text = element_text(family = "Roboto Slab"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
      ) +
      NULL

    ## only plot legend in case of more than 1 gene selected
    if(length(input$content) == 1){
      plot_complete  <-
        plot_complete +
        theme(legend.position = "none")
    }

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        labs(x = "") +
        theme(legend.position = "none")
    }

    return(plot_complete)
  }
  #error handling
  tryCatch(make_cellbar_raw(),
           error = function(x){make_bomb_plot()})
}

#figure legend
plot_cellbar_title <- "Dependency Bar Plot."
plot_cellbar_legend <- "Each bar shows the dependency scores of the queried genes in a cell line. Dependency scores less than -1 indicate a gene that is essential within a cell line. Dependency scores close to 0 mean no changes in fitness when the gene is knocked out. Dependency scores greater than 1 indicate gene knockouts lead to a gain in fitness."

#test
# make_cellbar(input = list(type = "gene", query = "ROCK1", content = c("ROCK1", "ROCK2")), mean = mean_virtual_achilles)


## DENSITY PLOT ----------------------------------------------------------------
#' Cell Dependencies Density Plot
#'
#' \code{make_cellbins} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a cellbins plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a cellbins plot. If an error is thrown, then will return a bomb plot.
#'
#' @export
#' @examples
#' make_cellbins(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_cellbins(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_cellbins(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_cellbins <- function(cellbins_data = achilles_long,
                          prism_data = prism_long,
                          expression_data = expression_names,
                          input = list(),
                          card = FALSE) {
  make_cellbins_raw <- function() {
    if(input$type == "gene") {
      xlab <- "Dependency Score (distribution)"
      aes_var <- rlang::sym("name")

      plot_data <-
        cellbins_data %>% #plot setup
        dplyr::filter(gene %in% input$content) %>%
        dplyr::left_join(expression_data, by = "X1") %>%
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
        prism_data %>% #plot setup
        dplyr::filter(name %in% input$content) %>%
        dplyr::left_join(expression_data, by = c("x1" = "X1")) %>%
        dplyr::select(-1) %>%
        dplyr::group_by(name) %>%
        dplyr::arrange(log2fc) %>%
        dplyr:: mutate(
          rank = 1:n(),
          med = median(log2fc, na.rm= TRUE)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::rename(dep_score = log2fc) #rename for graph

    } else if(input$type == "cell") {
      xlab <- "Dependency Scores (distribution)"
      aes_var <- rlang::sym("cell_line")

      plot_data <-
        cellbins_data %>% #plot setup
        dplyr::left_join(expression_data, by = "X1") %>%
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
      ggplot() +
      ## annotation range -1 to 1
      geom_rect(
        xmin = -1, xmax = 1,
        ymin = -Inf, ymax = Inf,
        fill = "grey95",
        show.legend = FALSE
      ) +
      ## indicator line y axis
      geom_linerange(
        aes(xmin = -Inf, xmax = med,
            y = fct_reorder(!!aes_var, -med),
            color = med < -1),
        linetype = "dotted",
        size = .2,
        show.legend = FALSE
      ) +
      ## density curves via {ggdist}
      stat_halfeye(aes(x = dep_score,
                       y = fct_reorder(!!aes_var, -med),
                       fill = stat(abs(x) > 1),
                       point_fill = after_scale(fill)),
                   .width = c(.025, .975),
                   color = "black",
                   shape = 21,
                   stroke = .7,
                   point_size = 2) +
      ## zero line
      geom_vline(
        xintercept = 0,
        color = "grey80",
        linetype = "dashed",
        show.legend = FALSE
      ) +
      ## titles
      labs(
        x = xlab,
        y = NULL,
        color = "Query",
        fill = "Query"
      ) +
      ## scales + legends
      scale_y_discrete(expand = c(.03, .03)) +
      scale_color_manual(values = c("grey70", colors)) +
      scale_fill_manual(values = c("grey70", colors)) +
      guides(
        color = guide_legend(size = 1, reverse = TRUE),
        fill = guide_legend(size = 1, reverse = TRUE)
      ) +
      ## theme changes
      theme_ddh() +
      theme(
        text = element_text(family = "Nunito Sans"),
        legend.position = "none",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(family = "Roboto Slab", size = 18),
        axis.text.x = element_text(size = 12, color = "grey30"),
        axis.title = element_text(size = 15)
      ) +
      NULL

    if(card) {
      plot_complete <-
        plot_complete +
        labs(x = "") + #too long of a xlabel
        NULL
    }
    return(plot_complete)
  }
  #error handling
  tryCatch(print(make_cellbins_raw()), #add print to catch ggplot.print errors
           error = function(x){make_bomb_plot()})
}

#figure legend
plot_cellbins_title <- "Computed Densities."
plot_cellbins_legend <- "Kernel density estimate of dependency scores. Dependency scores across all cell lines for queried genes, revealing overall influence of a gene on cellular fitness. The interval indicates the 95% quantile of the data, the dot indicates the median dependency score. The gray background highlights weak dependency values between -1 and 1."

#test
#make_cellbins(input = list(type = "gene", query = "ROCK1", content = c("ROCK1", "ROCK2")))
#make_cellbins(input = list(type = "compound", query = "aspirin", content = "aspiri"))
#make_cellbins(input = list(type = "cell", query = "HEPG2", content = "HEPG2"))

## LINEAGE LINERANGE PLOT ------------------------------------------------------
#' Dependency Lineage Plot
#'
#' \code{make_lineage} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a lineage plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a lineage plot. If an error is thrown, then will return a bomb plot.
#'
#' @export
#' @examples
#' make_lineage(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_lineage(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_lineage(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_lineage <- function(celldeps_data = achilles_long,
                         prism_data = prism_long,
                         expression_data = expression_names,
                         input = list(),
                         card = FALSE,
                         highlight = FALSE) {
  make_lineage_raw <- function() {
    if(input$type == "gene" ) {
      xlab <- "Dependency Score"
      title_var <- glue::glue('Cell lineage dependencies for {str_c(input$content, collapse = ", ")}')

      data_full <-
        celldeps_data %>% #plot setup
        dplyr::filter(gene %in% input$content) %>%
        dplyr::left_join(expression_data, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::mutate_at("lineage", function(str) {
          str <- str_replace_all(str, "\\_", " ")
          str <- str_to_title(str)
          return(str)
        }) %>%
        tidyr::drop_na(lineage) %>%
        tidyr::drop_na(dep_score) %>%
        dplyr::group_by(lineage) %>%
        dplyr::mutate(mean = mean(dep_score)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(lineage = fct_reorder(lineage, -mean))

      data_mean <- data_full %>%
        dplyr::group_by(lineage) %>%
        dplyr::summarize(dep_score = mean(dep_score))

    } else if(input$type == "compound") {
      xlab <- "Log2FC"
      title_var <- glue::glue('Cell lineage dependencies for {str_c(input$content, collapse = ", ")}')

      data_full <-
        prism_data %>% #plot setup
        dplyr::filter(name %in% input$content) %>%
        dplyr::left_join(expression_data, by = c("x1" = "X1")) %>%
        dplyr::select(-1) %>%
        dplyr::mutate_at("lineage", function(str) {
          str <- str_replace_all(str, "\\_", " ")
          str <- str_to_title(str)
          return(str)
        }) %>%
        tidyr::drop_na(lineage) %>%
        tidyr::drop_na(log2fc) %>%
        dplyr::group_by(lineage) %>%
        dplyr::mutate(mean = mean(log2fc)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(lineage = fct_reorder(lineage, -mean)) %>%
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
        group_by(lineage) %>%
        filter(n() > 1) %>%
        ungroup() %>%
        ggpubr::compare_means(dep_score ~ lineage, data = .,
                              ref.group = ".all.", method = "t.test") %>%
        filter(p.adj < 0.05)
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
      ggplot(aes(dep_score, lineage)) +
      ## annotation range -1 to 1
      # geom_rect(
      #   xmin = -1, xmax = 1,
      #   ymin = -Inf, ymax = Inf,
      #   fill = "grey95"
      # ) +
      ## zero line
      geom_vline(
        xintercept = 0,
        color = "grey80",
        linetype = "dashed"
      ) +
      ## indicator lines lineages
      geom_linerange(
        aes(xmin = -Inf, xmax = dep_score),
        color = "grey60",
        linetype = "dotted"
      ) +
      ## lineranges as "boxplots"
      stat_interval(
        data = data_full %>% {if(card)dplyr::filter(., lineage %in% most_negative) else .},
        orientation = "horizontal",
        .width = c(.05, .5, .95)
      ) +
      ## dot indicating mean
      geom_point(
        color = "black", fill = "white",
        shape = 21, stroke = .5,
        size = 1.8
      ) +
      ## scales + legends
      scale_x_continuous(
        sec.axis = dup_axis()
      ) +
      scale_color_ddh_d(
        palette = input$type,
        shuffle = TRUE, seed = 5L, ## to return "correctly" ordered, sequential colors
        labels = c("95%", "50%", "5%"),#of the data fall in these ranges
        name = ""
      ) +
      guides(color = guide_legend(reverse = TRUE)) +
      ## titles
      labs(
        x = xlab,
        y = NULL #,
        #title = title_var
      ) +
      {if(highlight)gghighlight(lineage %in% stats_data$group2, use_direct_label = FALSE)} + # toggle
      ## theme changes
      theme_ddh(grid = "none") +
      theme(
        legend.position = "top",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(family = "Roboto Slab"),
        axis.text.x = element_text(size = 12, color = "grey30"),
        axis.title = element_text(size = 15),
        axis.title.x.bottom = element_blank(),
        plot.title.position = "plot"
      ) +
      NULL

    if(card) {
      plot_complete <-
        plot_complete +
        scale_y_discrete(labels = scales::label_wrap(10)) + #to prevent long lines and squished plots
        scale_x_continuous(breaks = scales::breaks_extended(n = 3)) +
        theme(plot.title = element_blank())
    }
    return(plot_complete)
  }
  #error handling
  tryCatch(make_lineage_raw(),
           error = function(x){make_bomb_plot()})
}

#figure legend
plot_celllin_title <- "Cell Line Lineage Dependencies."
plot_celllin_legend <- "Each point shows the mean dependency score for the gene query within a given cell lineage. The intervals show the 5% quantiles centered on the median, the interquartile ranges, and the 95% quantiles. The gray background highlights weak dependency values between -1 and 1."

#make_lineage(input = list(type = "gene", content = c("ROCK2")))
#make_lineage(input = list(type = "gene", content = c("ROCK2")), card = TRUE)
#make_lineage(input = list(type = "gene", query = "ROCK1", content = c("ROCK1", "ROCK2")))
#make_lineage(input = list(type = "compound", query = "aspirin", content = "aspirin"))

## SUBLINE RANGE PLOT ---------------------------------------------------
#' Cell Dependency Sublineage Plot
#'
#' \code{make_sublineage} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a sublineage plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a sublineage plot. If an error is thrown, then will return a bomb plot.
#'
#' @export
#' @examples
#' make_sublineage(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_sublineage(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_sublineage(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_sublineage <- function(celldeps_data = achilles_long,
                            prism_data = prism_long,
                            expression_data = expression_names,
                            input = list(),
                            card = FALSE,
                            highlight = FALSE) {
  make_sublineage_raw <- function() {
    if(input$type == "gene") {
      xlab <- "Dependency Score"
      title_var <- glue::glue('Cell sub-lineage dependencies for {str_c(input$content, collapse = ", ")}')

      data_full <-
        celldeps_data %>% #plot setup
        dplyr::filter(gene %in% input$content) %>%
        dplyr::left_join(expression_data, by = "X1") %>%
        dplyr::select(-X1) %>%
        dplyr::mutate_at("lineage_subtype", function(str) {
          str <- str_replace_all(str, "\\_", " ")
          str <- if_else(str_detect(str, "^[:lower:]"), str_to_title(str), str)
          return(str)
        })  %>%
        tidyr::drop_na(lineage_subtype) %>%
        tidyr::drop_na(dep_score) %>%
        dplyr::group_by(lineage_subtype) %>%
        dplyr::mutate(mean = mean(dep_score)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(lineage_subtype = fct_reorder(lineage_subtype, -mean))

      data_mean <- data_full %>%
        dplyr::group_by(lineage_subtype) %>%
        dplyr::summarize(dep_score = mean(dep_score))

    } else if(input$type == "compound") {
      xlab <- "Log2FC"
      title_var <- glue::glue('Cell sub-lineage dependencies for {str_c(input$content, collapse = ", ")}')

      data_full <-
        prism_data %>% #plot setup
        dplyr::filter(name %in% input$content) %>%
        dplyr::left_join(expression_data, by = c("x1" = "X1")) %>%
        dplyr::select(-1) %>%
        dplyr::mutate_at("lineage_subtype", function(str) {
          str <- str_replace_all(str, "\\_", " ")
          str <- if_else(str_detect(str, "^[:lower:]"), str_to_title(str), str)
          return(str)
        })  %>%
        tidyr::drop_na(lineage_subtype) %>%
        tidyr::drop_na(log2fc) %>%
        dplyr::group_by(lineage_subtype) %>%
        dplyr::mutate(mean = mean(log2fc)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(lineage_subtype = fct_reorder(lineage_subtype, -mean)) %>%
        dplyr::rename(dep_score = log2fc) #rename for graph

      data_mean <- data_full %>%
        dplyr::group_by(lineage_subtype) %>%
        dplyr::summarize(dep_score = mean(dep_score))

    } else {
      stop("declare your type") }

    if(nrow(data_full) == 0){return(NULL)}

    if(highlight) {
      stats_data <- data_full %>%
        group_by(lineage_subtype) %>%
        filter(n() > 1) %>%
        ungroup() %>%
        ggpubr::compare_means(dep_score ~ lineage_subtype, data = .,
                              ref.group = ".all.", method = "t.test") %>%
        filter(p.adj < 0.05)
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
      ggplot(aes(dep_score, lineage_subtype)) +
      ## annotation range -1 to 1
      # geom_rect(
      #   xmin = -1, xmax = 1,
      #   ymin = -Inf, ymax = Inf,
      #   fill = "grey95"
      # ) +
      ## zero line
      geom_vline(
        xintercept = 0,
        color = "grey80",
        linetype = "dashed"
      ) +
      ## indicator lines sublineages
      geom_linerange(
        aes(xmin = -Inf, xmax = dep_score),
        color = "grey60",
        linetype = "dotted"
      ) +
      ## lineranges as "boxplots"
      stat_interval(
        data = data_full %>% {if(card)dplyr::filter(., lineage_subtype %in% most_negative) else .},
        orientation = "horizontal",
        .width = c(.05, .5, .95)
      ) +
      ## dot indicating mean
      geom_point(
        color = "black", fill = "white",
        shape = 21, stroke = .5,
        size = 1.8
      ) +
      ## scales + legends
      scale_x_continuous(
        sec.axis = dup_axis()
      ) +
      scale_color_ddh_d(
        palette = input$type,
        shuffle = TRUE, seed = 5L, ## to return "correctly" ordered, sequential colors
        labels = c("95%", "50%", "5%"), #of the data fall in these ranges
        name = ""
      ) +
      guides(color = guide_legend(reverse = TRUE)) +
      ## titles
      labs(
        x = xlab,
        y = NULL #,
        #title = title_var
      ) +
      {if(highlight)gghighlight(lineage_subtype %in% stats_data$group2, use_direct_label = FALSE)} + # toggle
      ## theme changes
      theme_ddh(grid = "none") +
      theme(
        legend.position = "top",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(family = "Roboto Slab"),
        axis.text.x = element_text(size = 12, color = "grey30"),
        axis.title = element_text(size = 15),
        axis.title.x.bottom = element_blank(),
        plot.title.position = "plot"
      ) +
      NULL

    if(card == TRUE) {
      plot_complete <-
        plot_complete +
        scale_y_discrete(labels = scales::label_wrap(10)) + #to prevent long lines and squished plots
        scale_x_continuous(breaks = scales::breaks_extended(n = 3)) +
        theme(plot.title = element_blank())
    }

    return(plot_complete)
  }
  #error handling
  tryCatch(make_sublineage_raw(),
           error = function(x){make_bomb_plot()})
}

#figure legend
plot_cellsublin_title <- "Cell Line Sub-Lineage Dependencies."
plot_cellsublin_legend <- "Each point shows the mean dependency score for the gene query within a given cell lineage. The intervals show the 5% quantiles centered on the median, the interquartile ranges, and the 95% quantiles. The gray background highlights weak dependency values between -1 and 1."

## CORRELATION PLOT FOR CELL DEPS--------------------------------------------------------
#' Co-essentiality Correlation Plot
#'
#' \code{make_correlation} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a correlation plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a correlation plot. If an error is thrown, then will return a bomb plot.
#'
#' @export
#' @examples
#' make_correlation(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_correlation(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_correlation(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_correlation <- function(table_data = achilles_cor_nest,
                             cell_line_data = achilles_cell_line_cor_nest,
                             prism_data = prism_cor_nest,
                             input = list(),
                             card = FALSE,
                             scale = NULL) { #no card option, but need this to prevent error
  make_correlation_raw <- function() {
    if(input$type == "gene") {
      mean <- mean_virtual_achilles
      upper_limit <- achilles_upper
      lower_limit <- achilles_lower
      var <- rlang::sym("fav_gene") #from https://rlang.r-lib.org/reference/quasiquotation.html
      label_var <- "Gene Rank"
      text_var <- "Gene"
      content_var <- glue::glue_collapse(input$content, sep = ", ")

      plot_data <-
        table_data %>%
        dplyr::filter(fav_gene %in% input$content) %>%
        tidyr::unnest(data) %>%
        dplyr::group_by(fav_gene) %>%
        dplyr::arrange(desc(r2)) %>%
        dplyr::mutate(
          rank = 1:n(),
          med = median(r2, na.rm= TRUE)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::rename(name = gene) #for plot name

    } else if(input$type == "compound") {
      mean <- mean_virtual_prism_cor
      upper_limit <- prism_cor_upper
      lower_limit <- prism_cor_lower
      var <- rlang::sym("fav_drug") #from https://rlang.r-lib.org/reference/quasiquotation.html
      label_var <- "Drug Rank"
      text_var <- "Compound"
      content_var <- glue::glue_collapse(input$content, sep = ", ")

      plot_data <-
        prism_data %>%
        dplyr::filter(fav_drug %in% input$content) %>%
        tidyr::unnest(data) %>%
        dplyr::group_by(fav_drug) %>%
        dplyr::arrange(desc(r2)) %>%
        dplyr::mutate(
          rank = 1:n(),
          med = median(r2, na.rm= TRUE)
        ) %>%
        dplyr::ungroup()

    }
    # else { #cell lines
    #   mean <- mean_virtual_achilles_cell_line
    #   upper_limit <- achilles_cell_line_upper
    #   lower_limit <- achilles_cell_line_lower
    #   var <- rlang::sym("fav_cell") #from https://rlang.r-lib.org/reference/quasiquotation.html
    #   label_var <- "Cell Rank"
    #   text_var <- "Cell"
    #   content_var <- glue::glue_collapse(input$content, sep = ", ")
    #
    #   plot_data <-
    #     cell_line_data %>%
    #     dplyr::filter(fav_cell %in% input$content) %>%
    #     tidyr::unnest(data) %>%
    #     dplyr::group_by(fav_cell) %>%
    #     dplyr::arrange(desc(r2)) %>%
    #     dplyr::mutate(
    #       rank = 1:n(),
    #       med = median(r2, na.rm= TRUE)
    #     ) %>%
    #     dplyr::ungroup() %>%
    #     dplyr::rename(name = cell) #for plot name
    # }

    if(!is.null(scale) & !card){
      plot_data <-
        plot_data %>%
        slice_sample(prop = scale)
    }

    if(card) {
      if(is.null(scale)) {
        scale <- 0.3
      }
      plot_data <- plot_data %>%
        group_by(name) %>%
        sample_n(scale*n()) %>% # scale is also used here
        ungroup()
    }

    plot_complete <-
      plot_data %>%
      ggplot() +
      ## sd square
      geom_hline(yintercept = upper_limit, color = "gray80") +
      geom_hline(yintercept = 0, color = "gray80", linetype = "dashed") +
      geom_hline(yintercept = lower_limit, color = "gray80") +
      annotate("text", x = Inf, y = 0.005, label = label_var, color = "gray40", hjust = 1.15 ,vjust = 0) +
      ## dot plot
      geom_point(aes(x = rank,
                     y = r2,
                     text = glue::glue('{text_var}: {name}'),
                     color = fct_reorder(!!var, med), #from https://rlang.r-lib.org/reference/quasiquotation.html
                     fill = fct_reorder(!!var, med) #from https://rlang.r-lib.org/reference/quasiquotation.html
      ),
      size = 1.1, stroke = .1, alpha = 0.4) +
      ## scales + legends
      scale_x_discrete(expand = expansion(mult = 0.02), na.translate = FALSE) +
      scale_color_ddh_d(palette = input$type) +
      scale_fill_ddh_d(palette = input$type) +
      guides(
        color = guide_legend(reverse = TRUE, override.aes = list(size = 5)),
        fill = guide_legend(reverse = TRUE, override.aes = list(size = 5))
      ) +
      ## labels
      labs(
        x = NULL,
        y = glue::glue("{text_var} correlations with {content_var}"),
        color = glue::glue("Query {text_var}"),
        fill = glue::glue("Query {text_var}")
      ) +
      ## theme changes
      theme_ddh() +
      theme(
        text = element_text(family = "Nunito Sans"),
        axis.text = element_text(family = "Roboto Slab"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
      ) +
      NULL

    ##change axis in case of more than 1 gene selected
    if(length(input$content) == 1){ #fix me
      plot_complete  <-
        plot_complete +
        guides(color = "none",
               fill = "none")
    }

    return(plot_complete)
  }
  #error handling
  tryCatch(make_correlation_raw(),
           error = function(x){make_bomb_plot()})
}

#figure legend
plot_genecorrelations_title <- "Correlation Curve."
#plot_genecorrelations_legend <- glue::glue("Each point shows the ranked correlation value ordered from high to low for each query. Correlation values outside the solid gray lines indicate the gene has a correlation value greater than {sd_threshold} standard deviations away from the mean.")
plot_genecorrelations_legend <- "Each point shows the ranked correlation value ordered from high to low for each query. Correlation values outside the solid gray lines indicate the gene has a correlation value greater than the mean."

#test
#make_correlation(input = list(type = "gene", query = "ROCK1", content = c("ROCK1")))
#make_correlation(input = list(type = "compound", query = "aspirin", content = "aspirin"))
#make_correlation(input = list(type = "cell", query = "HEPG2", content = "HEPG2"))

## EXPvDEP PLOT --------------------------------------------------------
#' Cell Expression v. Cell Dependency Plot
#'
#' \code{make_expdep} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a expdep plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a expdep plot. If an error is thrown, then will return a bomb plot.
#'
#' @export
#' @examples
#' make_expdep(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))
#' make_expdep(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)
#' \dontrun{
#' make_expdep(input = list(type = 'gene', content = 'ROCK1'))
#' }
make_expdep <- function(expression_data = expression_long,
                        celldeps_data = achilles_long,
                        expression_join = expression_names,
                        plot_se = TRUE,
                        input = list(),
                        card = FALSE) {
  make_expdep_raw <- function() {

    if (input$type == "gene") {
      exp_data <-
        expression_data %>% #plot setup
        dplyr::select(X1, gene, gene_expression) %>%
        dplyr::filter(gene %in% input$content)

      dep_data <-
        celldeps_data %>% #plot setup
        dplyr::filter(gene %in% input$content)

      combined_data <-
        exp_data %>%
        dplyr::inner_join(dep_data, by = c("X1", "gene")) %>%
        filter(!is.na(dep_score),
               !is.na(gene_expression)) %>%
        dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
        dplyr::mutate(med = median(dep_score, na.rm = TRUE)) %>%
        dplyr::left_join(expression_join, by = "X1") %>%
        dplyr::select(gene, gene_expression, dep_score, med, cell_line, lineage)
    }

    else if (input$type == "cell") {
      exp_data <-
        expression_data %>% #plot setup
        dplyr::left_join(expression_join, by = "X1") %>%
        dplyr::select(gene, cell_line, gene_expression) %>%
        dplyr::filter(cell_line %in% input$content)

      dep_data <-
        celldeps_data %>% #plot setup
        dplyr::left_join(expression_join, by = "X1") %>%
        dplyr::filter(cell_line %in% input$content)

      combined_data <-
        exp_data %>%
        dplyr::inner_join(dep_data, by = c("cell_line", "gene")) %>%
        filter(!is.na(dep_score),
               !is.na(gene_expression)) %>%
        dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>%
        dplyr::mutate(med = median(dep_score, na.rm = TRUE)) %>%
        dplyr::select(gene, gene_expression, dep_score, med, cell_line, lineage)
    }

    plot_complete <-
      combined_data %>%
      ggplot(aes(dep_score, gene_expression)) +
      ## gray background
      ## dot plot
      {if(input$type == "gene")geom_point(aes(color = fct_reorder(gene, med),
                                              fill = fct_reorder(gene, med)),
                                          size = 2, stroke = .1, alpha = 0.4)} +
      {if(input$type == "cell")geom_point(aes(color = fct_reorder(cell_line, med),
                                              fill = fct_reorder(cell_line, med)),
                                          size = 2, stroke = .1, alpha = 0.4)} +
      # smooth line
      {if(input$type == "gene")geom_smooth(aes(color = fct_reorder(gene, med),
                                               fill = fct_reorder(gene, med)),
                                           method = "lm",
                                           se = plot_se)} +
      {if(input$type == "cell")geom_smooth(aes(color = fct_reorder(cell_line, med),
                                               fill = fct_reorder(cell_line, med)),
                                           method = "lm",
                                           se = plot_se)} +
      # R coefs
      {if(card == FALSE & input$type == "gene")ggpubr::stat_cor(aes(color = fct_reorder(gene, med)),
                                                                digits = 3)} +
      {if(card == FALSE & input$type == "cell")ggpubr::stat_cor(aes(color = fct_reorder(cell_line, med)),
                                                                digits = 3)} +
      ## scales + legends
      scale_color_ddh_d(palette = input$type) +
      scale_fill_ddh_d(palette = input$type) +
      guides(
        color = guide_legend(reverse = TRUE, override.aes = list(size = 5)),
        fill = guide_legend(reverse = TRUE, override.aes = list(size = 5))
      ) +
      ## titles
      labs(
        x = "Dependency",
        y = "Expression",
        color = ifelse(input$type == "gene", "Query Gene", "Query Cell Line"),
        fill = ifelse(input$type == "gene", "Query Gene", "Query Cell Line")
      ) +
      ## theme changes
      theme_ddh() +
      theme(
        text = element_text(family = "Nunito Sans"),
        axis.text = element_text(family = "Roboto Slab")
      ) +
      NULL

    if(length(input$content) == 1){
      plot_complete  <-
        plot_complete +
        guides(color = "none",
               fill = "none") +
        labs(x = paste0(input$content, " Dependency"),
             y = paste0(input$content, " Expression"))
    }

    if(card == TRUE) {
      plot_complete <-
        plot_complete +
        labs(x = "", y = "") +
        theme(legend.position='none') +
        NULL
    }
    return(plot_complete)
  }
  #error handling
  tryCatch(make_expdep_raw(),
           error = function(x){make_bomb_plot()})
}

#test
#make_expdep(input = list(type = "gene", content = c("ROCK1")))
#make_expdep(input = list(type = "gene", content = c("ROCK1")), card = TRUE)
#make_expdep(input = list(type = "gene", content = c("ROCK1", "ROCK2")))
#make_expdep(input = list(type = "gene", content = c("ROCK1", "ROCK2")), card = TRUE)

#figure legend
plot_expdep_title <- "Gene Dependency versus Expression."
plot_expdep_legend <- paste0("Each point shows the dependency value compared to the expression value for gene within a given cell line. Gray area indicates dependency values that are between -1 and 1.")

#CELL --------------------------------------------------------------------
make_cell_image <- function(input = list()) {
  make_cell_image_raw <- function() {
    #if multiple, then pull single "rep" image; consider pulling >1 and using patchwork, eg.
    if(length(input$content > 1)) {
      fav_cell <- sample(input$content, 1) #input$gene_symbols[[1]]
    } else {
      fav_cell <- input$content
    }

    #fetch image path from dir, HARDCODED TO DATA DIR, NOT TEST DATA DIR
    file_name <- glue::glue('{fav_cell}_cell_image_plot.jpeg')
    image_path <- here::here(app_data_dir, "images", "cell", fav_cell, file_name)

    return(image_path)
  }
  #error handling
  tryCatch(make_cell_image_raw(),
           error = function(x){make_bomb_plot()})
}

#figure legend
plot_cellimage_title <- "Cell Images. "
plot_cellimage_legend <- "Cells image shown at low (top) and high (bottom) growth density."

#make_cell_image(input = list(content = "HEPG2"))

## CO-ESSENTIALITY CELL LINE PLOT --------------------------------------------------------
#' Cell Similarity Plot
#'
#' \code{make_cell_similarity_plot} returns an image of ...
#'
#' This is a plot function that takes a gene name and returns a cell similarity plot plot
#'
#' @param input Expecting a list containing type and content variable.
#' @param card A boolean that sets whether the plot should be scaled down to be a card
#' @return If no error, then returns a cell similarity plot If an error is thrown, then will return a bomb plot.
#'
#' @export
#' @examples
#' make_cell_similarity_plot(cell1 = "HEL", cell2 = "HEL9217")
#' make_cell_similarity_plot(cell1 = "HEL", cell2 = "KLE")
#' \dontrun{
#' make_cell_similarity_plot(cell1 = "HEL", cell2 = "HEPG2")
#' }
make_cell_similarity <- function(cell_sims_dep = cell_line_dep_sim,
                                 cell_sims_exp = cell_line_exp_sim,
                                 similarity = "dependency",
                                 input = list(),
                                 card = FALSE,
                                 scale = NULL) { #no card option, but need this to prevent error
  make_similarity_raw <- function() {

    if(similarity == "dependency") {
      cell_sims <- cell_sims_dep
    } else if(similarity == "expression") {
      cell_sims <- cell_sims_exp
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
      dplyr::arrange(desc(coef)) %>%
      dplyr::mutate(
        rank = 1:n(),
        med = median(coef, na.rm = TRUE)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::rename(name = cell1_name) #for plot name


    if(!is.null(scale) & !card){
      plot_data <-
        plot_data %>%
        slice_sample(prop = scale)
    }

    if(card) {
      if(is.null(scale)) {
        scale <- 0.3
      }
      plot_data <- plot_data %>%
        group_by(name) %>%
        sample_n(scale*n()) %>% # scale is also used here
        ungroup()
    }

    plot_complete <-
      plot_data %>%
      ggplot() +
      geom_hline(yintercept = 0, color = "gray80", linetype = "dashed") +
      annotate("text", x = Inf, y = 0.005, label = "Cell Rank", color = "gray40", hjust = 1.15 ,vjust = 0) +
      ## dot plot
      geom_point(aes(x = rank,
                     y = coef,
                     text = glue::glue('Cell: {name}'),
                     color = fct_reorder(name, med),
                     fill = fct_reorder(name, med)
      ),
      size = 1.1, stroke = .1, alpha = 0.4) +
      ## scales + legends
      scale_x_discrete(expand = expansion(mult = 0.02), na.translate = FALSE) +
      scale_color_ddh_d(palette = input$type) +
      scale_fill_ddh_d(palette = input$type) +
      guides(
        color = guide_legend(reverse = TRUE, override.aes = list(size = 5)),
        fill = guide_legend(reverse = TRUE, override.aes = list(size = 5))
      ) +
      ## labels
      labs(
        x = NULL,
        y = glue::glue("Cell coefficient estimates with {content_var}"),
        color = "Query Cell",
        fill = "Query Cell"
      ) +
      ## theme changes
      theme_ddh() +
      theme(
        text = element_text(family = "Nunito Sans"),
        axis.text = element_text(family = "Roboto Slab"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
      ) +
      NULL

    ##change axis in case of more than 1 gene selected
    if(length(input$content) == 1){ #fix me
      plot_complete  <-
        plot_complete +
        guides(color = "none",
               fill = "none")
    }

    return(plot_complete)
  }
  #error handling
  tryCatch(make_similarity_raw(),
           error = function(x){make_bomb_plot()})
}

#figure legend
make_cell_coessentiality_title <- "Cell Line Co-dependencies."
make_cell_coessentiality_legend <- glue::glue("Each point shows the ranked linear model coefficient estimate value ordered from high to low for each query.")

make_cell_coexpression_title <- "Cell Line Co-expressions."

# testing
# make_cell_similarity(input = list(type = "cell", query = "HEPG2", content = "HEPG2"))

# FUNCTIONAL PLOT ------------------------------------------
make_functional_cell <- function(pathway_data = pathways,
                                 expression_data = expression_long,
                                 expression_metadata = expression_meta,
                                 input = list(),
                                 num_pathways = 10,
                                 nwords = 5,
                                 num_genes = 2,
                                 remove_equivalent_pathways = FALSE,
                                 card = FALSE) {
  make_functional_cell_raw <- function() {

    plot_data <- pathway_data %>%
      unnest(data) %>%
      left_join(expression_data, by = "gene") %>%
      dplyr::select(pathway, go, gene, gene_expression, X1) %>%
      drop_na() %>%
      left_join(expression_metadata %>%
                  dplyr::select(X1, cell_line), by = "X1") %>%
      dplyr::select(-X1) %>%
      mutate(pathway_short = ifelse(str_count(pathway) >= nwords,
                                    paste0(gsub(paste0("^((\\w+\\W+){", nwords, "}\\w+).*$"), "\\1", pathway), " ..."),
                                    pathway)
      )

    med_cell <- plot_data %>%
      filter(cell_line %in% input$content) %>%
      group_by(go, cell_line) %>%
      mutate(genes_num = n()) %>%
      filter(genes_num >= num_genes) %>%
      ungroup() %>%
      group_by(cell_line, go) %>%
      summarise(med = median(gene_expression)) %>%
      ungroup() %>%
      filter(med != 0) %>% # decide if filter zero expressions
      left_join(plot_data %>%
                  dplyr::select(go, pathway_short) %>%
                  filter(!duplicated(go)),
                by = "go")

    if(remove_equivalent_pathways) {
      med_cell <- med_cell %>%
        filter(!duplicated(med))
    }

    if(length(input$content) == 1) {
      diff_cell <- plot_data %>%
        filter(!cell_line %in% input$content) %>%
        group_by(go) %>%
        mutate(genes_num = n()) %>%
        filter(genes_num > num_genes) %>%
        ungroup() %>%
        group_by(go) %>%
        summarise(med = median(gene_expression)) %>%
        ungroup() %>%
        filter(med != 0) %>% # decide if filter zero expressions
        bind_rows(med_cell) %>%
        group_by(go) %>%
        summarise(diff = abs(diff(med))) %>%
        ungroup() %>%
        dplyr::arrange(-diff) %>%
        dplyr::slice(1:num_pathways) %>%
        left_join(plot_data %>%
                    dplyr::select(go, pathway_short) %>%
                    filter(!duplicated(go)),
                  by = "go")
    } else {
      diff_cell <- med_cell %>%
        group_by(go) %>%
        summarise(diff = abs(diff(med))) %>%
        ungroup() %>%
        dplyr::arrange(-diff) %>%
        dplyr::slice(1:num_pathways) %>%
        left_join(plot_data %>%
                    dplyr::select(go, pathway_short) %>%
                    filter(!duplicated(go)),
                  by = "go")
    }

    med_cell <- med_cell %>%
      filter(go %in% diff_cell$go) %>%
      left_join(diff_cell %>%
                  dplyr::select(-pathway_short) %>%
                  filter(!duplicated(go)),
                by = "go")

    background <- plot_data %>%
      filter(go %in% diff_cell$go) %>%
      filter(!cell_line %in% input$content) %>%
      group_by(cell_line, go) %>%
      summarise(med = median(gene_expression)) %>%
      ungroup() %>%
      filter(med != 0) %>% # decide if filter zero expressions
      left_join(diff_cell, by = "go")

    plot_complete <- ggplot() +
      geom_jitter(data = background, aes(reorder(pathway_short, diff), med), color = "gray69", alpha = 0.5, width = 0.25) +
      geom_crossbar(data = med_cell,
                    aes(reorder(pathway_short, diff), med, ymin = med, ymax = med, color = cell_line)) +
      geom_point(data = med_cell, aes(reorder(pathway_short, diff), med, color = cell_line), size = 3) +
      labs(x = NULL,
           y = "Gene Expression",
           color = "Query Cell") +
      coord_flip() +
      scale_color_ddh_d(palette = input$type) +
      guides(
        color = guide_legend(reverse = TRUE, override.aes = list(size = 5))
      ) +
      ## theme changes
      theme_ddh() +
      theme(
        text = element_text(family = "Nunito Sans"),
        axis.text = element_text(family = "Roboto Slab"),
        legend.title = element_blank()
      ) +
      NULL

    if(length(input$content) == 1){
      plot_complete  <-
        plot_complete +
        ylab(paste0(str_c(input$content, collapse = ", "), " Gene Expression")) +
        theme(
          legend.position = "none"
        )
    } else {
      plot_complete <- plot_complete +
        theme(legend.title = element_blank())
    }

    if(card == TRUE){
      plot_complete <-
        plot_complete +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              legend.position = "none"
        )
    }

    return(plot_complete)
  }

  #error handling
  tryCatch(make_functional_cell_raw(),
           error = function(x){make_bomb_plot()})
}

functional_plot_title <- "Differential Pathway Expression Plot."
functional_plot_legend <- "The colored vertical bars indicate the pathway median expression for the queried cell line/s while the background grey points indicate the pathway median expression of all the other cell lines. If the query includes only one cell line, the difference between the median pathway expression of that cell line and the median pathway expression of all the other cell lines will be computed and the pathways with higher differences will appear first in the plot. Otherwise, the biggest difference between pathway medians of the queried cell lines will be used to rank the pathways in the plot. Those pathways with higher differences will appear first in the plot."

# testing
# make_functional_cell(input = list(type = "cell", content = c("HEL")))
# make_functional_cell(input = list(type = "cell", content = c("HEPG2", "HEL")))

# CELL METADATA PLOT ------------------------------------------
make_metadata_cell <- function(input = list(),
                               cell_line_similarity = "dependency",
                               metadata = "lineage",
                               bonferroni_cutoff = 0.05,
                               card = FALSE) {
  make_metadata_cell_raw <- function() {

    plot_data <- make_cell_sim_table(similarity = cell_line_similarity,
                                     bonferroni_cutoff = 1.1, # to include bonf == 1
                                     input = input) %>%
      bind_rows() %>%
      mutate(group = case_when(bonferroni > bonferroni_cutoff ~ "None",
                               bonferroni < bonferroni_cutoff & coef > 0 ~ "Similar",
                               bonferroni < bonferroni_cutoff & coef < 0 ~ "Dissimilar"),
             group = as_factor(group)
      ) %>%
      as_tibble()

    if(metadata == "lineage") {

      siglin <- plot_data %>%
        filter(bonferroni < bonferroni_cutoff) %>%
        pull(lineage) %>%
        unique()

      plot_data <- plot_data %>%
        filter(lineage %in% siglin)

      if(card & nrow(plot_data) > 6) {
        most_associated <-
          plot_data %>%
          arrange(pval) %>%
          filter(!duplicated(lineage)) %>%
          dplyr::slice(1:6) %>%
          dplyr::pull(lineage)

        plot_data <- plot_data %>%
          dplyr::filter(., lineage %in% most_associated)
      }

      plot_complete <- plot_data %>%
        filter(!is.na(lineage)) %>%
        group_by(lineage, group) %>%
        count() %>%
        ungroup() %>%
        ggplot(aes(x = n, y = lineage, fill = group, group = group))
    }
    else if (metadata == "sublineage") {

      sigsublin <- plot_data %>%
        filter(bonferroni < bonferroni_cutoff) %>%
        pull(lineage_subtype) %>%
        unique()

      plot_data <- plot_data %>%
        filter(lineage_subtype %in% sigsublin)

      plot_complete <- plot_data %>%
        filter(!is.na(lineage_subtype)) %>%
        group_by(lineage_subtype, group) %>%
        count() %>%
        ungroup() %>%
        ggplot(aes(x = n, y = lineage_subtype, fill = group, group = group))
    }

    plot_complete <- plot_complete +
      geom_col() +
      labs(x = "# Cell Lines",
           y = NULL) +
      # scale_x_continuous(labels = scales::percent_format()) +
      scale_fill_ddh_d(palette = input$type) +
      guides(
        color = guide_legend(reverse = TRUE, override.aes = list(size = 5))
      ) +
      ## theme changes
      theme_ddh() +
      theme(
        text = element_text(family = "Nunito Sans"),
        axis.text = element_text(family = "Roboto Slab"),
        legend.title = element_blank(),
        legend.position = "top"
      ) +
      NULL

    if(card) {
      plot_complete <-
        plot_complete +
        scale_y_discrete(labels = scales::label_wrap(10)) + #to prevent long lines and squished plots
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              plot.title = element_blank(),
              legend.position = "none"
        )
    }

    return(plot_complete)
  }

  #error handling
  tryCatch(make_metadata_cell_raw(),
           error = function(x){make_bomb_plot()})
}

cell_metadata_plot_title <- "Lineage Similarity Plot."
cell_metadata_plot_legend <- "Similar and dissimilar lineages (and sublineages) associated with the queried cell line/s."

# testing
# make_metadata_cell(input = list(type = "cell", content = c("HEL")))
# make_metadata_cell(input = list(type = "cell", content = c("HEPG2", "HEL")))

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
#' @export
#' @examples
#' make_molecule_structure(input = list(content = "aspirin"))
#' make_molecule_structure(input = list(content = "ibuprofen"))
#' make_molecule_structure(input = list(content = 2244))
#' \dontrun{
#' make_molecule_structure(input = list(content = "aspirin"))
#' }
make_molecule_structure <- function(input = list(),
                                    #test = FALSE,
                                    model_width = 800,
                                    model_height = 800,
                                    card = FALSE) {
  make_molecule_structure_raw <- function(){
    #Customize the material with toon shading
    shiny_toon_material = rayvertex::material_list(type="toon_phong",
                                                   toon_levels=,
                                                   toon_outline_width=0.1)

    molecule <-
      raymolecule::get_molecule(input$content) %>%
      raymolecule::generate_full_scene(pathtrace=FALSE,
                                       material_vertex = shiny_toon_material) %>%
      raymolecule::render_model(width=model_width,
                                height=model_height,
                                background="white")

    if(card == TRUE){
      molecule <-
        molecule +
        #annotate with image magick?
        labs(x = "")#, title = "Gene Information", caption = "more ...")
    }

    return(molecule)
  }
  #error handling
  tryCatch(make_molecule_structure_raw(),
           error = function(x){make_bomb_plot()})
}

#make_molecule_structure(input = list(content = "aspirin"))
#make_molecule_structure(content = 2244)
#make_molecule_structure(content = "ibuprofen")
#make_molecule_structure(content = "ibuprofe")

