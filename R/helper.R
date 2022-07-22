
#' Function to create an empty
#'
#' @return A ggplot2 object.
#'
#' @examples
#' make_empty_plot()
make_empty_plot <- function() {
  ggplot() +
    labs(title = "Nothing to see here. Try again.")
}

#' Function to create a bomb plot
#'
#' @return A ggplot2 object.
#'
#' @examples
#' make_bomb_plot()
make_bomb_plot <- function(){
  #fill background
  background <- data.frame(
    x = rep(1:72),
    y = rep(1:72, each = 72)
  )

  #set circle w/randomness (strictly controlled)
  x_center <- round(sample(x = 28:32,
                           size = 1), digits = 0) #30
  y_center <- round(sample(x = 28:32,
                           size = 1), digits = 0) #30
  radius <- round(sample(x = 18:22,
                         size = 1), digits = 0) #20

  #drop out those points outside of circle
  # https://www.redblobgames.com/grids/circle-drawing/
  circle <- data.frame(
    x = rep((x_center-radius):(radius+x_center)),
    y = rep((y_center-radius):(radius+y_center), each = length(rep((x_center-radius):(radius+x_center))))
  )

  circle <-
    circle %>%
    mutate(distance_to_center = sqrt((x-x_center)^2 + (y-y_center)^2),
           include = if_else(distance_to_center < radius, TRUE, FALSE))

  #add bomb top
  #target 1/3 of bomb width; target 1/6 for height
  bomb_top_width <- (radius*2)/4
  bomb_top_height <- bomb_top_width/3
  bomb_top <- data.frame(
    x = rep((x_center-(bomb_top_width)/2):(x_center+(bomb_top_width)/2)),
    y = rep((y_center+radius):(y_center+radius+bomb_top_height), each = length(rep((x_center-(bomb_top_width)/2):(x_center+(bomb_top_width)/2))))
  )

  #add white streak
  streak_arc <- radius*.75
  white_streak <- data.frame(
    y = c(rep((y_center+1):(y_center+streak_arc)),rep((y_center-1):(y_center-streak_arc)))
  )

  white_streak <-
    white_streak %>%
    mutate(x = round(30 - (sqrt((streak_arc^2) - (y_center-y)^2)), digits = 0))

  #drop 2/3
  samples_to_keep <- nrow(white_streak)/3

  #build fuse
  cube <- sample(1:2, size = 1) #2
  lin <- sample(1:10, size = 1) #3
  whole <- 1 #2

  fuse_fun <- function(x,
                       cube_var=2,
                       lin_var=3,
                       whole_var=2){
    y = cube_var*(x^3) - lin_var*(x) + whole_var #- 3*(x)^2
    return(y)
  }

  fuse <-
    data.frame(
      x = seq(from = -1,
              to = 2.5,
              length.out = radius))

  fuse <-
    fuse %>%
    mutate(y = fuse_fun(x,
                        cube_var = cube,
                        lin_var = lin,
                        whole_var = whole))

  # fuse <-
  #   fuse %>%
  #   mutate(y1 = fuse_fun(x),
  #          y2 = fuse_fun(x, whole_var = 20),
  #          y3 = fuse_fun(x, whole_var = 2),
  #          y4 = fuse_fun(x, whole_var = 0.1))
  #
  # ggplot() +
  #   geom_line(data = fuse, aes (x, y1), color = "blue") +
  #   geom_line(data = fuse, aes (x, y2), color = "red") +
  #   geom_line(data = fuse, aes (x, y3), color = "green") +
  #   geom_line(data = fuse, aes (x, y4), color = "pink")

  fuse <-
    fuse %>%
    dplyr::mutate(x_final = rep(x_center:(x_center+radius-1)),
                  y_final = round((y_center + radius + bomb_top_height + y), digits = 0)) %>%
    dplyr::filter(x_final < 66, y_final < 66)

  #add first point to fuse
  fuse_start <- #add 2 vertical blocks to start fuse
    data.frame(
      x_final = x_center,
      y_final = round(y_center + radius + bomb_top_height + 1, digits = 0))

  #find last point in fuse
  x_fuse_bind <-
    fuse %>%
    dplyr::slice_max(x_final) %>%
    dplyr::pull(x_final)

  y_fuse_bind <-
    fuse %>%
    dplyr::slice_max(x_final) %>%
    dplyr::pull(y_final)

  fuse_end <- #add 3 horizontal blocks to fuse
    data.frame(
      x_final = c(x_fuse_bind, x_fuse_bind + 1, x_fuse_bind + 2),
      y_final = c(y_fuse_bind, y_fuse_bind, y_fuse_bind))

  #bind rows
  fuse <-
    fuse %>%
    dplyr::bind_rows(fuse_start) %>%
    dplyr::bind_rows(fuse_end)

  #spark
  #find last point in fuse
  x_fuse <-
    fuse %>%
    dplyr::slice_max(x_final) %>%
    dplyr::pull(x_final)

  y_fuse <-
    fuse %>%
    dplyr::slice_max(x_final) %>%
    dplyr::pull(y_final)

  #find points around a circle of set size
  spark_radius <- round(radius/3, digits=0)

  spark <- data.frame(
    x = rep((x_fuse-spark_radius):(spark_radius+x_fuse)),
    y = rep((y_fuse-spark_radius):(spark_radius+y_fuse), each = length(rep((x_fuse-spark_radius):(spark_radius+x_fuse))))
  )

  spark <-
    spark %>%
    dplyr::mutate(distance_to_center = sqrt((x-x_fuse)^2 + (y-y_fuse)^2),
                  include = dplyr::if_else(distance_to_center < spark_radius, TRUE, FALSE)) %>%
    dplyr::filter(x < 72,
                  y < 72)

  spark_diameter <- spark_radius * 2
  spark_thirds <- round(spark_diameter/3, digits = 0)

  #select random x positions
  x_left_pos = round(sample(x = (x_fuse-spark_radius):(x_fuse-spark_radius+spark_thirds),
                            size = 1), digits=0)
  x_middle_pos = round(sample(x = (x_fuse-spark_thirds/2):(x_fuse+spark_thirds/2),
                              size = 1), digits=0)
  x_right_pos = round(sample(x = (x_fuse+spark_radius-spark_thirds):(x_fuse+spark_radius),
                             size = 1), digits=0)
  x_left_neg = round(sample(x = (x_fuse-spark_radius):(x_fuse-spark_radius+spark_thirds),
                            size = 1), digits=0)
  x_middle_neg = round(sample(x = (x_fuse-spark_thirds/2):(x_fuse+spark_thirds/2),
                              size = 1), digits=0)
  x_right_neg = round(sample(x = (x_fuse+spark_radius-spark_thirds):(x_fuse+spark_radius),
                             size = 1), digits=0)

  #calculate y positions
  sparky <- function(x,
                     x_origin = x_fuse,
                     y_origin = y_fuse,
                     pos = TRUE){
    if(pos == TRUE){
      y_offset = sqrt((spark_radius^2 - (abs(x_origin - x))^2))
      y = round(y_origin + y_offset, digits = 0)
    } else {
      y_offset = sqrt((spark_radius^2 - (abs(x_origin - x))^2))
      y = round(y_origin - y_offset, digits = 0)
    }
    return(y)
  }

  y_left_pos = round(sparky(x_left_pos), digits=0)
  y_middle_pos = round(sparky(x_middle_pos), digits=0)
  y_right_pos = round(sparky(x_right_pos), digits=0)
  y_left_neg = round(sparky(x_left_neg, pos = FALSE), digits=0)
  y_middle_neg = round(sparky(x_middle_neg, pos = FALSE), digits=0)
  y_right_neg = round(sparky(x_right_neg, pos = FALSE), digits=0)

  #for testing
  spark_points <- data.frame(
    x = c(x_left_pos, x_middle_pos, x_right_pos, x_left_neg, x_middle_neg, x_right_neg),
    y = c(y_left_pos, y_middle_pos, y_right_pos, y_left_neg, y_middle_neg, y_right_neg)
  )

  spark_line_maker <- function(x_origin = x_fuse,
                               y_origin = y_fuse,
                               x_var,
                               y_var){
    #add slope finding fun
    slope_finder <- function(x1 = x_origin,
                             x2 = x_var,
                             y1 = y_origin,
                             y2 = y_var){
      m = (y2-y1)/(x2-x1)
      return(m)
    }
    #get slope
    slope <- slope_finder()

    #make first df
    spark_frame <-
      data.frame(
        y = rep(y_origin:y_var)
      )
    final_frame <-
      spark_frame %>%
      dplyr::mutate(x = round((y-y_origin)/slope, digits = 0)+x_origin)
    return(final_frame)
  }

  sparkline1 <- spark_line_maker(x_var = x_left_pos, y_var = y_left_pos)
  sparkline2 <- spark_line_maker(x_var = x_middle_pos, y_var = y_middle_pos)
  sparkline3 <- spark_line_maker(x_var = x_right_pos, y_var = y_right_pos)
  sparkline4 <- spark_line_maker(x_var = x_left_neg, y_var = y_left_neg)
  sparkline5 <- spark_line_maker(x_var = x_middle_neg, y_var = y_middle_neg)
  sparkline6 <- spark_line_maker(x_var = x_right_neg, y_var = y_right_neg)

  spark_lines <-
    sparkline1 %>%
    bind_rows(sparkline2) %>%
    bind_rows(sparkline3) %>%
    bind_rows(sparkline4) %>%
    bind_rows(sparkline5) %>%
    bind_rows(sparkline6)

  #drop some
  spark_lines <-
    spark_lines %>%
    mutate(distance_to_center = sqrt((x-x_fuse)^2 + (y-y_fuse)^2),
           include = if_else(distance_to_center < spark_radius/3, FALSE, TRUE))

  #build plot with layers
  plot_complete <-
    ggplot() +
    geom_tile(data = background, aes(x, y), fill = "white", color = "white") +
    geom_tile(data = circle %>% dplyr::filter(include == TRUE), aes(x, y), fill = "black", color = "black") +
    geom_tile(data = bomb_top, aes(x, y), fill = "black", color = "black") +
    geom_tile(data = white_streak %>% sample_n(size = samples_to_keep), aes(x, y), fill = "white", color = "white") +
    geom_tile(data = fuse, aes(x_final, y_final), fill = "black", color = "black") +
    #geom_tile(data = spark %>% dplyr::filter(include == TRUE), aes(x, y), fill = "blue", alpha = 0.5) +
    #geom_tile(data = spark_points, aes(x, y), fill = "blue") +
    geom_tile(data = spark_lines %>% dplyr::filter(include == TRUE), aes(x, y), fill = "black", color = "black") +
    coord_cartesian(xlim = c(0,72), ylim = c(0,72)) +
    theme_void() +
    NULL
  return(plot_complete)
}

#' Function to obtain plot sizes
#'
#' @param function_name Function name.
#'
#' @return A list with plot size
#'
#' @examples
#' plot_size_finder("make_ideogram")
plot_size_finder <- function(function_name){ #this function sets the output size to pre-render a plot
  #switch statement for majority of plots, ignore for cards only, and COMMENT out intentionally excluded
  type <-
    switch(function_name,
           make_ideogram = "portrait",
           make_proteinsize = "landscape",
           #make_sequence = "square",
           #make_protein_domain_plot = "landscape",
           make_radial = "square",
           #make_umap_plot = "landscape",
           #make_cluster_enrich_plot = "square",
           make_structure = "portrait",
           make_pubmed =  "landscape",
           make_cellanatogram =  "landscape",
           make_female_anatogram = "portrait",
           make_male_anatogram = "portrait",
           make_tissue = "landscape",
           make_celldeps = "landscape",
           make_cellbins = "landscape",
           make_lineage = "portrait",
           make_sublineage = "portrait",
           make_correlation = "landscape",
           make_cellexpression = "landscape",
           make_cellgeneprotein = "landscape",
           make_expdep = "landscape",
           stop("no such plot")

    )

  if(type == "portrait"){
    #Size: 1080 x 1920 px
    #Aspect Ratio: 9:16
    # plot_size = list(plot_width = 1080,
    #                  plot_height = 1920)
    # plot_size = list(plot_width = 480,
    #                  plot_height = 720)
    plot_size = list(plot_width = 300,
                     plot_height = 450)
  } else if(type == "landscape") {
    #Size: 1920 x 1080 px
    #Aspect ratio: 16:9
    # plot_size = list(plot_width = 1920,
    #                  plot_height = 1080)
    # plot_size = list(plot_width = 720,
    #                  plot_height = 480)
    plot_size = list(plot_width = 750,
                     plot_height = 400)
  } else if(type == "square"){
    #Aspect ratio: 1:1
    plot_size = list(plot_width = 500,
                     plot_height = 500)
  } else {
    #custom?
    plot_size = list(plot_width = 720,
                     plot_height = 720)
  }
  return(plot_size)
}

#' Function to create roxygen structures
#'
#' @param fun_name Function name.
#' @param type Function type.
#'
#' @return A roxygen structure.
#'
#' @examples
#' make_roxygen("make_ideogram")
make_roxygen <- function(fun_name,
                         type = "graph"){
  name <- stringr::str_remove_all(fun_name, pattern = "make_") %>% str_replace_all(pattern = "_", replacement = " ")
  title <- stringr::str_to_title(glue::glue("{name} {type}"))
  roxygen_c <- c(
    glue::glue("#' {title}"),
    "#'",
    glue::glue("#' \\code{{{fun_name}}} returns an image of ..."),
    "#'",
    glue::glue("#' This is a {type} function that takes a gene name and returns a {name} {type}"),
    "#'",
    "#' @param input Expecting a list containing type and content variable.",
    if(type == "plot"){"#' @param card A boolean that sets whether the plot should be scaled down to be a card"},
    glue::glue("#' @return If no error, then returns a {name} {type}. If an error is thrown, then will return an empty {type}."),
    "#'",
    "#' @export",
    "#' @examples",
    glue::glue("#' {fun_name}(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'))"),
    if(type == "plot"){glue::glue("#' {fun_name}(input = list(type = 'gene', query = 'ROCK1', content = 'ROCK1'), card = TRUE)")},
    "#' \\dontrun{",
    glue::glue("#' {fun_name}(input = list(type = 'gene', content = 'ROCK1'))"),
    "#' }"
    )
  roxygen <- unlist(roxygen_c)
  output_file <- here::here("tmp.Rmd")
  writeLines(roxygen, con = output_file)
}

