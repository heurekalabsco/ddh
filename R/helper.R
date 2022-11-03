#' Function to download all DDH .Rds files
#'
#' @param app_data_dir Data directory path.
#' @param object_name Optional object name to get a single file; default gets all files in bucket
#' @param test Boolean that will download test data instead of full data
#' @param overwrite Boolean that deletes the app_data_dir thereby removing files, before remaking dir and downloading all files
#' @param privateMode Boolean indicating if private data is required.
#'
#' @importFrom magrittr %>%
#'
#' @export
download_ddh_data <- function(app_data_dir,
                              object_name = NULL,
                              test = FALSE,
                              private = TRUE,
                              overwrite = FALSE
){
  #delete dir to overwrite
  if(overwrite == TRUE) {
    if(is.null(object_name)){all_objects <- ""} else {all_objects <- object_name}
    path <- stringr::str_glue("{app_data_dir}/{all_objects}")
    path %>% purrr::walk(unlink, recursive = TRUE) #handles 1, many, or NULL file_name
  }
  #make dir
  if(!dir.exists(app_data_dir)){
    dir.create(app_data_dir)
  }

  #get_data function
  get_aws_data <- function(bucket, #public, private, test
                           object_name){
    bucket_var <-
      switch(bucket,
             public = "AWS_DATA_BUCKET_ID",
             private = "AWS_DATA_BUCKET_ID_PRIVATE",
             # raw = "AWS_DATA_BUCKET_ID_RAW",
             test = "AWS_DATA_BUCKET_ID_TEST")

    s3 <- paws::s3()
    data_objects <-
      s3$list_objects(Bucket = Sys.getenv(bucket_var)) %>%
      purrr::pluck("Contents")

    message(glue::glue('{length(data_objects)} objects in the {bucket} bucket'))

    #filter list for single object
    if(!is.null(object_name)){
      key_name <- glue::glue('_data/{object_name}')
      data_objects <-
        data_objects %>% #take full list
        purrr::keep(purrr::map_lgl(.x = 1:length(data_objects),
                                   ~ data_objects[[.x]][["Key"]] %in% key_name)) #pass map_lgl to keep to filter names to keep
      message(glue::glue('filtered to keep {length(data_objects)}'))
    }

    #check for no objects
    if(length(data_objects) == 0){
      return(message(glue::glue("No objects found. Check your object name.")))
    }

    for (i in 1:length(data_objects)) {
      key_name <- data_objects[[i]][["Key"]]
      file_name <- stringr::str_remove(data_objects[[i]][["Key"]], "_data/") #remove subfolder structure

      #check if files exists
      if(file.exists(glue::glue("{app_data_dir}/{file_name}"))){
        message(glue::glue("file already exists: {file_name}"))
      } else {
        #if not, then download
        s3$download_file(Bucket = Sys.getenv(bucket_var),
                         Key = as.character(key_name),
                         Filename = as.character(glue::glue("{app_data_dir}/{file_name}")))
        message(glue::glue("file downloaded: {file_name}"))
      }
    }
  }

  #test data
  if(test == TRUE){
    get_aws_data(bucket = "test",
                 object_name = object_name)
    return(message("test data download complete"))
  }
  #get private data
  if(private == TRUE){
    get_aws_data(bucket = "private",
                 object_name = object_name)
    #no return
  }

  #get data
  get_aws_data(bucket = "public",
               object_name = object_name)
  return(message("data download complete"))
}

#' Function to Load All DDH Data Including .Rds Files and Colors
#'
#' @param app_data_dir Data directory path.
#' @param file_name Optional file name to get a single file; default loads all files
#' @param privateMode Boolean indicating if private data is required.
#'
#' @export
load_ddh_data <- function(app_data_dir,
                          object_name = NULL,
                          feather = TRUE,
                          db = FALSE) {

  # Load colors
  load_ddh_colors()
  message("loaded colors")

  if(feather == TRUE) {
    # Load .RDS files
    load_ddh_feather(app_data_dir,
                     object_name)
    message("loaded Rds files")
  }

  if(!is.null(object_name)){ #stop here
    return(message("done"))
  }

  if(db) {
    #load DB cons
    load_ddh_db()
    message("loaded db connections")
  }

  message("finished loading")
}

#' Function to load all DDH files
#'
#' @param app_data_dir Data directory path.
#' @param object_name Optional object name to load a single file; default loads all files
#'
#' @importFrom magrittr %>%
#'
#' @export
load_ddh_feather <- function(app_data_dir,
                             object_name = NULL) {
  if(is.null(object_name)){ #all objects
    all_objects <- fs::dir_ls(path = app_data_dir)
  } else if(object_name %in% c("gene", "cell", "compound")) { #object "group"
    object_regex <- stringr::str_c(object_name, "universal", sep = "|")
    all_objects <-
      fs::dir_ls(path = app_data_dir) %>%
      stringr::str_subset(pattern = object_regex)
  } else { #single object
    all_objects <- glue::glue('{app_data_dir}/{object_name}')
  }

  #file loader constructor
  load_feather_object <- function(filename) {
    object <- basename(filename) #removes filepath
    object <- sub("_test", "", object)
    assign(object, arrow::read_feather(filename), envir = .GlobalEnv)

    message(glue::glue("loaded {object}"))
  }

  #walk through to load all files
  all_objects %>%
    purrr::walk(load_feather_object)

  #print done
  message("finished loading feathers")
}

#' Function to load all DDH db connections
#'
#' @param object_name Optional object name to load a single db connection; default loads all connections
#'
#' @importFrom magrittr %>%
#'
#' @export
load_ddh_db <- function(object_name = NULL) {
  #fun for db connection
  #load_db_connection <- function(db_name){}

  #manually (?) list all dbs to connect to

  #filter list for those in file_name
  #db_cons <-
  # if(is.null(file_name)){
  #   db_cons <-  #manually (?) list all dbs to connect to
  # } else {
  #   db_cons <- file_name
  # }

  #db_cons %>% purrr::walk(load_db_connection)

  #placeholder for db connections
  message("placeholder for db connections")
}

#' Function to load DDH colors
#'
#' @export
load_ddh_colors <- function() {
  ## MAIN COLORS
  ##2EC09C  ## cyan
  ##BE34EF  ## violet
  ##E06B12  ## orange
  ##004AAB  ## blue
  ##F0CE44  ## yellow
  ##1785A4  ## blend cyan + blue

  load_colors <- function() {
    ## Color sets  for genes
    color_set_gene <- generate_colors("#2EC09C")  ## cyan
    #CC is the hex alpha conversion for 80%, so the next line adds it; used in graph
    color_set_gene_alpha <- purrr::map_chr(color_set_gene, ~ glue::glue_collapse(c(.x, "CC"), sep = ""))

    ## Palette function for genes
    color_pal_gene <- grDevices::colorRampPalette(color_set_gene)

    ## Color sets for proteins
    color_set_protein <- generate_colors("#004AAB")  ## blue
    color_set_protein_alpha <- purrr::map_chr(color_set_protein, ~ glue::glue_collapse(c(.x, "CC"), sep = ""))

    ## Palette function for proteins
    color_pal_protein <- grDevices::colorRampPalette(color_set_protein)

    ## Color sets for proteins
    color_set_geneprotein <- generate_colors("#1785A4")  ## blue
    color_set_geneprotein_alpha <- purrr::map_chr(color_set_geneprotein, ~ glue::glue_collapse(c(.x, "CC"), sep = ""))

    ## Palette function for proteins
    color_pal_geneprotein <- grDevices::colorRampPalette(color_set_geneprotein)

    ## Color sets for cells
    color_set_cell <- generate_colors("#BE34EF")  ## violet
    color_set_cell_alpha <- purrr::map_chr(color_set_cell, ~ glue::glue_collapse(c(.x, "CC"), sep = ""))

    ## Palette function for cells
    color_pal_cell <- grDevices::colorRampPalette(color_set_cell)

    ## Color sets for compounds
    color_set_compound <- generate_colors("#E06B12")  ## orange
    color_set_compound_alpha <- purrr::map_chr(color_set_compound, ~ glue::glue_collapse(c(.x, "CC"), sep = ""))

    ## Palette function for compounds
    color_pal_compound <- grDevices::colorRampPalette(color_set_compound)

    return(list(color_set_gene=color_set_gene,
                color_set_gene_alpha=color_set_gene_alpha,
                color_pal_gene=color_pal_gene,
                color_set_protein=color_set_protein,
                color_set_protein_alpha=color_set_protein_alpha,
                color_pal_protein=color_pal_protein,
                color_set_geneprotein=color_set_geneprotein,
                color_set_geneprotein_alpha=color_set_geneprotein_alpha,
                color_pal_geneprotein=color_pal_geneprotein,
                color_set_cell=color_set_cell,
                color_set_cell_alpha=color_set_cell_alpha,
                color_pal_cell=color_pal_cell,
                color_set_compound=color_set_compound,
                color_set_compound_alpha=color_set_compound_alpha,
                color_pal_compound=color_pal_compound)
    )
  }

  ddh_colors <- load_colors()
  list2env(ddh_colors, .GlobalEnv)

}

#' Function to create an empty table
#'
#' @return A data.frame.
#' @export
#'
#' @examples
#' make_empty_table()
make_empty_table <- function() {
  message("No data available")
}

#' Function to create an empty plot
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' make_empty_plot()
make_empty_plot <- function() {
  ggplot2::ggplot() +
    ggplot2::labs(title = "Nothing to see here. Try again.")
}

#' Function to create a bomb plot
#'
#' @importFrom magrittr %>%
#'
#' @return A ggplot2 object.
#' @export
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
    dplyr::mutate(distance_to_center = sqrt((x-x_center)^2 + (y-y_center)^2),
                  include = dplyr::if_else(distance_to_center < radius, TRUE, FALSE))

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
    dplyr::mutate(x = round(30 - (sqrt((streak_arc^2) - (y_center-y)^2)), digits = 0))

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
    dplyr::mutate(y = fuse_fun(x,
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
    dplyr::bind_rows(sparkline2) %>%
    dplyr::bind_rows(sparkline3) %>%
    dplyr::bind_rows(sparkline4) %>%
    dplyr::bind_rows(sparkline5) %>%
    dplyr::bind_rows(sparkline6)

  #drop some
  spark_lines <-
    spark_lines %>%
    dplyr::mutate(distance_to_center = sqrt((x-x_fuse)^2 + (y-y_fuse)^2),
                  include = dplyr::if_else(distance_to_center < spark_radius/3, FALSE, TRUE))

  #build plot with layers
  plot_complete <-
    ggplot2::ggplot() +
    ggplot2::geom_tile(data = background, ggplot2::aes(x, y), fill = "white", color = "white") +
    ggplot2::geom_tile(data = circle %>% dplyr::filter(include == TRUE), ggplot2::aes(x, y), fill = "black", color = "black") +
    ggplot2::geom_tile(data = bomb_top, ggplot2::aes(x, y), fill = "black", color = "black") +
    ggplot2::geom_tile(data = white_streak %>% dplyr::sample_n(size = samples_to_keep), ggplot2::aes(x, y), fill = "white", color = "white") +
    ggplot2::geom_tile(data = fuse, ggplot2::aes(x_final, y_final), fill = "black", color = "black") +
    #ggplot2::geom_tile(data = spark %>% dplyr::filter(include == TRUE), ggplot2::aes(x, y), fill = "blue", alpha = 0.5) +
    #ggplot2::geom_tile(data = spark_points, ggplot2::aes(x, y), fill = "blue") +
    ggplot2::geom_tile(data = spark_lines %>% dplyr::filter(include == TRUE), ggplot2::aes(x, y), fill = "black", color = "black") +
    ggplot2::coord_cartesian(xlim = c(0,72), ylim = c(0,72)) +
    ggplot2::theme_void() +
    NULL
  return(plot_complete)
}

#' Function to obtain plot sizes
#'
#' @param function_name Function name.
#'
#' @importFrom magrittr %>%
#'
#' @return A list with plot size
#'
#' @export
#' @examples
#' plot_size_finder("make_ideogram")
plot_size_finder <- function(function_name){ #this function sets the output size to pre-render a plot
  #switch statement for majority of plots, ignore for cards only, and COMMENT out intentionally excluded
  type <-
    switch(function_name,
           make_ideogram = "portrait",
           make_proteinsize = "short",
           #make_sequence is handeled in generate_sequences()
           #make_protein_domain_plot = "landscape",
           make_radial = "square",
           #make_umap_plot = "landscape",
           #make_cluster_enrich_plot = "square",
           make_pubmed =  "landscape",
           make_cellanatogram =  "landscape",
           make_female_anatogram = "portrait",
           make_male_anatogram = "portrait",
           make_tissue = "custom_tissue",
           make_celldeps = "landscape",
           make_cellbins = "landscape",
           make_lineage = "tall",
           make_sublineage = "taller",
           make_correlation = "landscape",
           make_cellexpression = "short",
           make_cellgeneprotein = "short",
           make_expdep = "landscape"
           #stop("no such plot") #invisibly returns NULL for first conditional
    )
  if(is.null(type)){
    #setting to NA let's ggplot decide on scales
    plot_size = list(plot_width = NA,
                     plot_height = NA)
  } else if(type == "portrait"){
    #Size: 1080 x 1920 px
    #Aspect Ratio: 9:16
    plot_size = list(plot_width = 1080,
                     plot_height = 1920)
  } else if(type == "landscape") {
    #Size: 1920 x 1080 px
    #Aspect ratio: 16:9
    plot_size = list(plot_width = 1920,
                     plot_height = 1080)
  } else if(type == "square"){
    #Aspect ratio: 1:1
    plot_size = list(plot_width = 1920,
                     plot_height = 1920)
  } else if(type == "custom_tissue"){
    #Aspect ratio: 1:1
    plot_size = list(plot_width = 1920,
                     plot_height = 3840) #1920*2
  } else if(type == "short"){
    plot_size = list(plot_width = 1920,
                     plot_height = 720) #1080/1.5
  } else if(type == "tall"){
    #lineage
    plot_size = list(plot_width = 2160, #1080*2
                     plot_height = 2496) #1920*1.3
  } else if(type == "taller"){
    #sublineage
    plot_size = list(plot_width = 2872, #1080*2.66
                     plot_height = 4473) #1920*2.33
  } else {
    #setting to NA let's ggplot decide on scales
    plot_size = list(plot_width = NA,
                     plot_height = NA)
  }
  return(plot_size)
}

#' Function to create roxygen structures
#'
#' @param fun_name Function name.
#' @param type Function type.
#'
#' @importFrom magrittr %>%
#'
#' @return A roxygen structure.
#'
#' @export
#' @examples
#' make_roxygen("make_ideogram")
make_roxygen <- function(fun_name,
                         type = "graph"){
  name <- stringr::str_remove_all(fun_name, pattern = "make_") %>%
    stringr::str_replace_all(pattern = "_", replacement = " ")
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

#' Empty Graph Graph
#'
#' \code{make_empty_graph} returns an image of ...
#'
#' This is a graph function that takes a gene name and returns a empty graph graph
#'
#' @param input Expecting a list containing type and content variable.
#' @return If no error, then returns a empty graph graph. If an error is thrown, then will return an empty graph.
#'
#' @importFrom magrittr %>%
#' @import visNetwork
#'
#' @export
#' @examples
#' make_empty_graph()
#' \dontrun{
#' make_empty_graph()
#' }
make_empty_graph <- function(type = "gene") {
  if(type == "gene") {
    queryColor <- color_set_gene_alpha[2]
  } else if (type == "cell") {
    queryColor <- color_set_cell_alpha[2]
  }else if (type == "compound") {
    queryColor <- color_set_compound_alpha[2]
  } else {
    stop("declare your type")
  }
  #set params, copied from below
  borderColor <- "rgba(204, 204, 204, 0.8)" #(gray80), formerly white 255, 255, 255
  displayHeight = '90vh'
  displayWidth = '100%'
  #make empty nodes table
  nodes_empty = dplyr::tibble(id = 0,
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

#' Extract Documentation Sections
#'
#' Extract parts of a function help documentation.
#'
#' @param fun Function name.
#' @param section Help section to extract.
#'
#' @export
#' @examples
#' help_extract("make_radial", package = ddh, section = "Description")
#' help_extract("make_radial", package = ddh, section = "Examples")
help_extract <- function(fun,
                         section = "Description",
                         ...) {
  x <- capture.output(tools:::Rd2txt(utils:::.getHelpFile(help(fun, ...)),
                                     options = list(sectionIndent = 0)))
  B <- grep("^_", x)
  x <- gsub("_\b", "", x, fixed = TRUE)
  X <- rep(FALSE, length(x))
  X[B] <- 1
  out <- split(x, cumsum(X))
  out <- out[[which(sapply(out, function(x)
    grepl(section, x[1], fixed = TRUE)))]][-c(1, 2)]
  while(TRUE) {
    out <- out[-length(out)]
    if (out[length(out)] != "") { break }
  }

  if(section == "Description") {
    out <- paste0(out, collapse = " ")
  }

  return(out)
}

#' Extract Title Function
#'
#' Extract function title from documentation.
#'
#' @param fun Function name.
#'
#' @export
#' @examples
#' title_extract("make_radial", package = ddh)
#' title_extract("make_umap_plot", package = ddh)
title_extract <- function(fun,
                          ...) {
  x <- capture.output(tools:::Rd2txt(utils:::.getHelpFile(help(fun, ...)),
                                     options = list(sectionIndent = 0)))
  x <- gsub("_\b", "", x, fixed = TRUE)
  title <- x[1]

  return(title)
}

#' Make Legend
#'
#' @param fun Function name.
#' @param html Boolean indicating HTML format or not. Default is FALSE.
#'
#' @export
#' @examples
#' make_legend(fun = "make_radial")
make_legend <- function(fun,
                        html = FALSE,
                        ...) {
  description <- help_extract(fun, package = ddh, section = "Description")
  title <- title_extract(fun, package = ddh)

  if(html) {
    legend <- paste0("<b>", title, ".</b> ", description)
  } else {
    legend <- paste0(title, ". ", description)
  }

  return(legend)
}

#DATA GENERATION----
#' Fix names
#'
#' @param wrong_name A gene name that is not its official symbol
#' @param summary_df The gene_summary dataframe with all gene symbols
#'
#' @importFrom magrittr %>%
#'
#' @export
fix_names <- function(wrong_name,
                      summary_df = gene_summary) {
  var <- stringr::str_which(summary_df$aka, paste0("(?<![:alnum:])", wrong_name, "(?![:alnum:]|\\-)")) #finds index
  df <- summary_df[var,]
  right_name <- df$approved_symbol
  if (length(var) == 1) {
    return(right_name)
  } else {
    return(wrong_name)
  }
  #fixes 251, leaves 11
}

#' Clean colnames
#'
#' @param dataset A dataset to fix gene names
#' @param summary_df The gene_summary dataframe with all gene symbols
#'
#' @importFrom magrittr %>%
#'
#' @export
clean_colnames <- function(dataset,
                           summary_df = gene_summary) {
  for (name in names(dataset)) {
    if (name %in% summary_df$approved_symbol == FALSE) {
      fixed_name <- fix_names(name, summary_df = summary_df)
      if (fixed_name %in% names(dataset) == FALSE) {
        names(dataset)[names(dataset) == name] <- fixed_name
      }
    }
  }
  return(dataset)
}

#' Get Stats
#'
#' @param data_universal_stats_summary Dataframe of stats summary generated in ddh_data.
#' @param data_set A character indicating data set from which the stats were generated, typically one of achilles, expression_gene, expression_protein, or prism name.
#' @param var A character indicating the variable to extract from the stats summary dataframe, typically one of threshold, sd, mean, upper, pr lower.
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' get_stats(data_set = "achilles", var = "sd")
get_stats <- function(data_universal_stats_summary = universal_stats_summary,
                      data_set,
                      var){
  stat <-
    data_universal_stats_summary %>%
    dplyr::filter(id == data_set) %>%
    dplyr::pull(var)
}

#CARD HELPERS----
#' Load single card
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' load_image(input = list(type = "gene", content = c("ROCK3")), fun_name = "make_female_anatogram")
#' load_image(input = list(type = "gene", content = c("ROCK1")), fun_name = "make_female_anatogram", card = TRUE)
#' load_image(input = list(type = "gene", content = c("ROCK1", "ROCK2")), fun_name = "make_female_anatogram")
#' load_image(input = list(type = "compound", content = c("aspirin")), fun_name = "make_celldeps")
#' load_image(input = list(type = "compound", content = c("aspirin")), fun_name = "make_molecule_structure")
load_image <- function(input = list(),
                       fun_name,
                       card = FALSE) { #type is either card or plot
  load_image_raw <- function(){
    #build the input for the raw fun() & call function
    name <- stringr::str_c(input$content, collapse="-") #intended to fail with multigene query to return NULL
    fun <- stringr::str_remove(fun_name, "make_")
    if(card == TRUE){image_type = "card"} else {image_type = "plot"}

    file_name <- glue::glue('{name}/{name}_{fun}_{image_type}.jpg')

    #check to see if file exists
    #check if exists
    url <- glue::glue("https://{Sys.getenv('AWS_IMAGES_BUCKET_ID')}.s3.amazonaws.com/{file_name}")
    status <- httr::GET(url) %>% httr::status_code()

    if(status == 200){
      return(url)
    } else {
      return(NULL)
    }
  }
  #error handling
  tryCatch(load_image_raw(),
           error = function(e){
             message(e)
             NULL
           })
}

#' Load PDB file
#'
#' @param input Expecting a list containing type and content variable.
#'
#' @return Path to a url containing a PDB file.
#'
#' @export
#' @examples
#' load_pdb(input = list(content = c("ROCK1")))
#' load_pdb(input = list(content = c("ROCK3")))
#' load_pdb(input = list(content = c("ROCK1", "ROCK2")))
load_pdb <- function(input = list()){
  load_pdb_raw <- function(){
    if(length(input$content > 1)){
      gene_symbol <- input$content[1]
    } else {
      gene_symbol <- input$content
    }

    #check if exists
    url <- glue::glue("https://{Sys.getenv('AWS_PROTEINS_BUCKET_ID')}.s3.amazonaws.com/{gene_symbol}.pdb")
    status <- httr::GET(url) %>% httr::status_code()

    if(status == 200){
      return(url)
    } else {
      return(NULL)
    }
  }
  #error handling
  tryCatch(load_pdb_raw(),
           error = function(e){
             message(e)
             NULL
           })
}

#' Format Path Part
#'
#' @export
format_path_part <- function(key) {
  # formats part of an image path replacing invalid characters slashes with "_"
  gsub("[/]", "_", key)
}

#QUARTO HELPER----
#' Make Quarto
#'
#' @importFrom magrittr %>%
#'
#' @export
#' @examples
#' make_quarto("Gene Information")
make_quarto <- function(title){
  file_name <- janitor::make_clean_names(title) %>%
    stringr::str_replace_all(pattern = "_", replacement = "-")
  quarto_c <- c(
    glue::glue("# {title} {{#sec-{file_name}}}"),
    "  ",
    "The purpose of this page is to provide...  ",
    "The data comes from...  ")

  quarto_c <- unlist(quarto_c)
  file_path <- here::here("code", "methods")
  output_file <- glue::glue("{file_path}/{file_name}.qmd")
  writeLines(quarto_c, con = output_file)

  return(glue::glue("{file_name}.qmd"))
}

#' Good File Namer
#'
#' @param input Expecting a list that contains a content object
#'
#' @export
#'
#' @examples
#' good_file_namer(input = list(content = c("ROCK1", "ROCK2")))
good_file_namer <- function(input = list){
  if(length(input$content) == 1){
    good_file_name <- input$content
  } else if(length(input$content) > 1){
    good_file_name <- paste0("custom_query_", paste0(stringr::str_extract(input$content, pattern = "^[:alnum:]+"),collapse = "_"))
  }
  return(good_file_name)
}

#INTERNAL LINKS----
#' Internal link
#'
#' @importFrom magrittr %>%
#'
#' @export
internal_link <- function(query, linkout_img=FALSE) {
  hrefr <- function(x, linkout_img) {
    paste0('<a href="?show=gene&query=',
           x,
           '" target="_blank">',
           dplyr::if_else(linkout_img==TRUE, '<img src="link out_25.png", width="10", height="10">', x),
           '</a>')
  }
  if (stringr::str_detect(query, ",") == TRUE){ #code for list
    string_vec <- unlist(stringr::str_split(query, ", "))
    string_vec_link <- purrr::map_chr(string_vec, hrefr, linkout_img = FALSE)
    query_link <- stringr::str_c(string_vec_link, collapse = ", ")

  } else if (stringr::str_detect(query, "^[:digit:]{7}") == TRUE) { #regex for GO, used in shiny_tables browsePathwaysPanelServer
    query_link <- paste0('<a href="?show=pathway&go=',
                         query,
                         '">',
                         query,
                         '</a>')
  } else {
    query_link <- hrefr(query, linkout_img)
  }
  return(query_link)
}

#' Internal link cell
#'
#' @importFrom magrittr %>%
#'
#' @export
internal_link_cell <- function(query) {
  query_link <-
    paste0('<a href="?show=cell&cell_line=',
           query,
           '" target="_blank">',
           query,
           '</a>')

  return(query_link)
}

#' PubMed Linkr
#'
#' @importFrom magrittr %>%
#'
#' @export
pubmed_linkr <- function(query, number_only = FALSE) { #add for make_pubmed_table()
  if (stringr::str_detect(query, "PubMed:[:digit:]{1,9}|PubMed=[:digit:]{1,9}") == TRUE) {
    num <- stringr::str_extract(query, "[:digit:]{1,9}")
    link <- paste0("https://pubmed.ncbi.nlm.nih.gov/", num)
    href_link <- paste0('<a href="',
                        link,
                        '" target="_blank">', #open in new page
                        num,
                        '</a>')
    query_link <- stringr::str_replace(query, "PubMed:[:digit:]{1,9}", href_link)
    return(query_link)
  } else if (stringr::str_detect(query, "[:digit:]{1,9}") == TRUE && number_only == TRUE) {
    num <- query
    link <- paste0("https://pubmed.ncbi.nlm.nih.gov/", num)
    href_link <- paste0('<a href="',
                        link,
                        '" target="_blank">', #open in new page
                        num,
                        '</a>')
    query_link <- stringr::str_replace(query, "[:digit:]{1,9}", href_link)
    return(query_link)
  } else {
    return(query)
  }
}

#' PMC Linkr
#'
#' @importFrom magrittr %>%
#'
#' @export
pmc_linkr <- function(query) { #add for make_pubmed_table()
  link <- paste0("https://www.ncbi.nlm.nih.gov/pmc/articles/", query)
  href_link <- paste0('<a href="',
                      link,
                      '" target="_blank">', #open in new page
                      query,
                      '</a>')
  query_link <- stringr::str_replace(query, "PMC[:digit:]{1,9}", href_link)
  return(query_link)
}

#' Uniprot Linkr
#'
#' @importFrom magrittr %>%
#'
#' @export
uniprot_linkr <- function(query) {
  if (stringr::str_detect(query, "UniProtKB:[:alnum:]{1,6}") == TRUE) {
    num <- stringr::str_extract(query, "UniProtKB:[:alnum:]{1,6}") %>%
      stringr::str_split(pattern = ":")  %>%
      purrr::pluck(1, 2)
    link <- paste0("https://www.uniprot.org/uniprot/", num)
    href_link <- paste0('<a href="',
                        link,
                        '" target="_blank">', #open in new page
                        num,
                        '</a>')
    query_link <- stringr::str_replace(query, "UniProtKB:[:alnum:]{1,6}", href_link)
    return(query_link)
  } else {
    return(query)
  }
}

#' Uniprot Linkr 2
#'
#' @importFrom magrittr %>%
#'
#' @export
uniprot_linkr2 <- function(query) {
  link <- paste0("https://www.uniprot.org/uniprot/", query)
  query_link <- paste0('<a href="',
                       link,
                       '" target="_blank">',
                       query,
                       '</a>')
  return(query_link)
}

#' PDB Linkr
#'
#' @importFrom magrittr %>%
#'
#' @export
pdb_linkr <- function(query) {
  link <- paste0("https://www.rcsb.org/structure/", query)
  query_link <- paste0('<a href="',
                       link,
                       '" target="_blank">',
                       query,
                       '</a>')
  return(query_link)
}

#' Eco Zapr
#'
#' @importFrom magrittr %>%
#'
#' @export
eco_zapr <- function(query) {
  if (stringr::str_detect(query, "ECO:[:digit:]{1,7}") == TRUE) {
    query_zap <- stringr::str_remove_all(query, "ECO:[:digit:]{1,7}\\|")
    query_zap <- stringr::str_remove_all(query_zap, "ECO:[:digit:]{1,7}")
    return(query_zap)
  }
  else {
    return(query)
  }
}

#' Bracketr
#'
#' @importFrom magrittr %>%
#'
#' @export
bracketr <- function(query) {
  if (stringr::str_detect(query, "\\{") == TRUE) {
    query_plus <- stringr::str_replace(query, "\\{", "\\{Curated links: ")
    return(query_plus)
  }
  else {
    return(query)
  }
}

#' Gene Linkr
#'
#' @importFrom magrittr %>%
#'
#' @export
gene_linkr <- function(summary_table = gene_summary, query) { #best to use with lit_linkr, which takes a char_vec
  query_clean <- stringr::str_extract(query, "[^[:punct:]]+") #anything but punct, one or more
  if ((query_clean %in% summary_table$approved_symbol) == TRUE) {
    query_link <- paste0('<a href="?show=gene&query=',
                         query_clean,
                         '" target="_blank">',
                         query,
                         '</a>')
    return(query_link)
    # NOT YET WORKING; AKA SEARCH EITHER RETURNS NO MATCHES B/C QUERY =! MULTIPLE GENES IN STR OR QUERY == TOO MANY GENES IN UNLISTED STR
    # } else if (str_detect(unlist(str_split(gene_summary$aka, ", ")), query_clean) == TRUE) { #search unlisted AKA
    #   aka_num <- str_which(summary_table$aka, query_clean)
    #   query_aka <- summary_table[aka_num, 1]
    #   query_link <- paste0('<a href="?show=gene&query=',
    #                        query_aka,
    #                        '" target="_blank">',
    #                        query,
    #                        '</a>')
    #   return(query_link)
  } else {
    return(query)
  }
}

#' Lit Linkr
#'
#' @importFrom magrittr %>%
#'
#' @export
lit_linkr <- function(summary_table = gene_summary,
                      lit_string) { #split into char_vec, map function, glue back together, and then treat as htmlOutput in shiny_text
  lit_string <- stringr::str_replace_all(lit_string, "PubMed ", "PubMed:") #if space, then str_split breaks
  lit_string <- stringr::str_replace_all(lit_string, "PubMed=", "PubMed:") #for CC
  lit_vec <- unlist(stringr::str_split(lit_string, pattern = " "))
  lit_vec <- purrr::map_chr(lit_vec, eco_zapr) #remove ECO
  lit_links <-
    lit_vec %>%
    purrr::map_chr(pubmed_linkr) %>% #make pubmed links
    purrr::map_chr(uniprot_linkr) %>%  #make uniprot links
    purrr::map_chr(bracketr) %>% #annotate bracket
    purrr::map_chr(gene_linkr, summary_table = gene_summary) #add some internal links
  lit_string_link <- stringr::str_c(lit_links, collapse = " ")
  return(lit_string_link)
}

#' Drug Linkr
#'
#'
#' @importFrom magrittr %>%
#'
#' @export
drug_linkr <- function(query) {
  if ((query %in% prism_names$name) == TRUE) {
    query_link <- paste0('<a href="?show=compound&query=',
                         query,
                         '" target="_blank">',
                         query,
                         '</a>')
    return(query_link)
  } else {
    return(query)
  }
}

#' MOA Linkr
#'
#'
#' @importFrom magrittr %>%
#'
#' @export
moa_linkr <- function(query) {
  if ((query %in% prism_names$moa) == TRUE) {
    query_link <- paste0('<a href="?show=moa&query=',
                         stringr::str_replace_all(query, " ", "%20"),
                         '" target="_blank">',
                         stringr::str_to_title(query),
                         '</a>')
    return(query_link)
  } else {
    return(query)
  }
}

#' Metabolite Linkr
#'
#'
#' @importFrom magrittr %>%
#'
#' @export
metabolite_linkr <- function(query) {
  if ((query %in% hmdb_names$name) == TRUE) {
    query_link <- paste0('<a href="?show=compound&query=', #fix me
                         query,
                         '" target="_blank">',
                         query,
                         '</a>')
    return(query_link)
  } else {
    return(query)
  }
}

#' Cell Linkr
#'
#' Make cell, lineage, list linkr with type var
#'
#' @importFrom magrittr %>%
#'
#' @export
cell_linkr <- function(query, type) {
  type_url <- switch (type,
                      cell = "?show=cell&cell_line=",
                      lineage = "?show=lineage&query=",
                      lineage_subtype = "?show=lineage_subtype&query="
  )
  if (query %in% expression_names$cell_line | query %in% expression_names$lineage | query %in% expression_names$lineage_subtype) {
    query_no_space <- stringr::str_replace_all(query, " ", "%20")
    query_title <- query #stringr::str_to_title(query)
    query_link <- glue::glue('<a href="{type_url}{query_no_space}" target="_blank">{query_title}</a>')
    return(query_link)
  } else {
    return(query)
  }
}

###make example
#' Make Example
#'
#' Function to create an HTML string of examples to use on the index page and the start here methods doc
#'
#' @param privateMode Boolean indicating if private data is pointed to
#' @return An HTML string.
#'
#' @export
#' @examples
#' make_example()
#' \dontrun{
#' make_example()
#' }
make_example <- function(privateMode = TRUE){
  #need to dynamically switch depending on privateMode
  tld <- dplyr::if_else(privateMode == TRUE, "com", "org")
  url <- glue::glue('http://www.datadrivenhypothesis.{tld}/')

  #need full urls so works in methods, instead of reference url
  examples <- glue::glue('<h4>Search for genes</h5>
              <ul>
                <li>A single gene, such as <a href="{url}?show=gene&query=TP53">TP53</a> or <a href="{url}?show=gene&query=BRCA1">BRCA1</a></li>
                <li>A pathway name, such as <a href="{url}?show=search&query=cholesterol">cholesterol</a>, which will lead you to <a href="{url}?show=pathway&query=0006695">Cholesterol Biosynthetic Process</a></li>
                <li>The Gene Ontology biological process identifier, such as <a href="{url}?show=search&query=1901989">1901989</a>, which will find <a href="{url}?show=pathway&query=1901989">Pathway: Positive Regulation Of Cell Cycle Phase Transition (GO:1901989)</a></li>
                <li>A custom list of genes (separated by commas), such as <a href="{url}?show=search&query=BRCA1,%20BRCA2">BRCA1, BRCA2</a>, which will search <a href="{url}?show=gene_list&query=BRCA1,BRCA2">a custom gene list</a></li>
              </ul>')
  return(examples)
}

