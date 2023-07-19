# LOAD DDH COLORS-----------------------------------------------------
#'
#' Loads and defines colors & palettes for DDH
#'
#' @return Returns a list of color sets and their alpha versions along with palette functions for genes
#' proteins, cells, and compounds, which can then be loaded into DDH
#'
#' @examples
#' data_colors <- load_colors()
#'
#' @author Matthew Hirschey & Pol Castellano
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

#' GENERATE COLORS----------------------------------------------
#'
#' Generates color palettes using a specified main color
#'
#' @param hex Hex code for main color as string
#' @param discrete If true, categorical data is reordered to place main color first
#'
#' @return A set, or palette, of three colors
#'
#' @examples
#' generate_colors("#2EC09C")
#' generate_colors("#BE34EF")
#'
#' @author Matthew Hirschey & Pol Castellano
#'
#' @export
generate_colors <- function(hex, discrete = FALSE) {
  color_set <- c(
    colorspace::desaturate(colorspace::lighten(hex, .6), .2),
    hex,
    colorspace::darken(hex, .6, space = "HLS")
  )

  if (discrete == TRUE) color_set <- color_set[c(2, 1, 3)]

  return(color_set)
}


#' COLOR RAMP PALETTE SHUFFLE ----------------------------------
#'
#' Generates a shuffled color palette based on `colorRampPalette`
#'
#' @param colors Colors to interpolate; must be a valid argument to col2rgb().
#' @param seed If NULL, set to a randomized integer. Otherwise, takes on given integer value
#'
#' @return A function with an argument vector containing values between 0 and 1 that are
#'         mapped to a numeric matrix of RGB color values with one row per color
#'         and 3 or 4 columns
#'
#' @examples
#' colorRampPalette_shuffle(colors = c("#DEBBF4", "#BE34EF", "#AC7B84"))
#' colorRampPalette_shuffle(colors = c("#DEBBF4", "#BE34EF", "#AC7B84"), seed = 123L)
#'
#' @author Matthew Hirschey & Pol Castellano
#'
#' @export
colorRampPalette_shuffle <- function (colors, seed = NULL, ...)
{
  if (is.null(seed)) seed <- sample.int(999999999, 1)
  if (!is.integer(seed)) stop('seed should be an integer value.')
  set.seed(seed)

  ramp <- colorRamp(colors, ...)
  function(n) {
    x <- ramp(seq.int(0, 1, length.out = n))
    rand <- sample(nrow(x))
    x <- x[rand, ]
    if (ncol(x) == 4L)
      rgb(x[, 1L], x[, 2L], x[, 3L], x[, 4L], maxColorValue = 255)
    else rgb(x[, 1L], x[, 2L], x[, 3L], maxColorValue = 255)
  }
}

#' DDH COLORS --------------------------------------------------
#'
#' Extracts DDH colors as hex codes
#'
#' @param ... Categorical character names of colors - some combination of "gene", "protein",
#' "gene_protein", "cell", and "compound"
#'
#' @return If no arguments provided, DDH's full color set is returned. If parameters provided,
#' the corresponding categories' color subsets are returned.
#'
#' @examples
#' ddh_colors()
#' ddh_colors("gene")
#'
#' @author Matthew Hirschey & Pol Castellano
#'
#' @export
ddh_colors <- function(...) {
  ddh_cols <- c(
    `gene`         = "#2EC09C",
    `protein`      = "#004AAB",
    `gene_protein` = "#1785A4",
    `cell`         = "#BE34EF",
    `compound`     = "#024223"
  )

  cols <- c(...)

  if (is.null(cols))
    return (ddh_cols)

  ddh_cols[cols]
}


#' DDH PAL C
#'
#' Returns a function that can easily generate and interpolate a continuous DDH color palette
#'
#' @param palette Character name of palette in ddh_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments to pass to [grDevices::colorRampPalette()]
#'
#' @return A color palette function created using `colorRampPalette`
#'
#' @examples
#' ddh_pal_c()(10)
#' ddh_pal_c("protein")(3)
#' ddh_pal_c("compound", reverse = TRUE)(5)
#'
#' @author Matthew Hirschey & Pol Castellano
#'
#' @export
ddh_pal_c <- function(palette = "gene", reverse = FALSE, ...) {

  palette <- stringr::str_to_lower(palette)

  if (!palette %in% c("gene", "protein", "gene_protein", "cell", "compound")) stop('palette should be one of "gene", "protein", "gene_protein", "cell", or "compound".')

  ddh_palettes <- list(
    `gene`         = unname(generate_colors(ddh_colors("gene"))),
    `protein`      = unname(generate_colors(ddh_colors("protein"))),
    `gene_protein` = unname(generate_colors(ddh_colors("gene_protein"))),
    `cell`         = unname(generate_colors(ddh_colors("cell"))),
    `compound`     = unname(generate_colors(ddh_colors("compound")))
  )

  pal <- ddh_palettes[[palette]]

  if (reverse) pal <- rev(pal)

  grDevices::colorRampPalette(pal, ...)
}


#' DDH PAL D
#'
#' Returns a function that can interpolate a discrete DDH color palette
#'
#' @param palette Character name of palette in ddh_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param shuffle Boolean indicating whether the colors should be shuffled
#' @param seed Boolean setting a seed to ensure consistent randomization
#'
#' @return A color palette function created using `colorRampPalette`
#'
#' @examples
#' ddh_pal_d()(10)
#' ddh_pal_d("protein")(3)
#' ddh_pal_d("compound", reverse = TRUE, shuffle = FALSE)(5)
#'
#' @author Matthew Hirschey & Pol Castellano
#'
#' @export

ddh_pal_d <- function(palette = "gene", reverse = FALSE, shuffle = FALSE, seed = NULL) {

  palette <- stringr::str_to_lower(palette)

  if (is.null(seed)) seed <- sample.int(999999999, 1)

  if (!palette %in% c("gene", "protein", "gene_protein", "cell", "compound")) stop('palette should be one of "gene", "protein", "gene_protein", "cell", or "compound".')
  if (!is.logical(reverse)) stop('reverse should be logical.')
  if (!is.logical(shuffle)) stop('shuffle should be logical.')
  if (!is.integer(seed)) stop('seed should be an integer value.')

  ddh_palettes <- list(
    `gene`         = unname(generate_colors(ddh_colors("gene"), discrete = TRUE)),
    `protein`      = unname(generate_colors(ddh_colors("protein"), discrete = TRUE)),
    `gene_protein` = unname(generate_colors(ddh_colors("gene_protein"), discrete = TRUE)),
    `cell`         = unname(generate_colors(ddh_colors("cell"), discrete = TRUE)),
    `compound`     = unname(generate_colors(ddh_colors("compound"), discrete = TRUE))
  )

  pal <- ddh_palettes[[palette]]

  if (reverse) pal <- rev(pal)

  if (shuffle == TRUE) {
    colorRampPalette_shuffle(pal, seed = seed)
  } else {
    grDevices::colorRampPalette(pal)
  }
}


#' SCALE COLOR DDH C
#'
#' Color scale constructor for continuous DDH color palettes
#'
#' @param palette Character name of palette in ddh_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE
#'            or FALSE
#'
#' @return A color scale for continuous color palettes for various categories
#'
#' @examples
#' library(ggplot2)
#' ggplot(mpg, aes(hwy, cty, color = displ)) + geom_point(size = 4) +
#' scale_color_ddh_c()
#' ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Sepal.Width)) +
#' geom_point(size = 4) + scale_color_ddh_c("cOMPOUND", reverse = TRUE)
#'
#' @author Matthew Hirschey & Pol Castellano
#'
#' @export
scale_color_ddh_c <- function(palette = "gene", reverse = FALSE, ...) {

  palette <- stringr::str_to_lower(palette)

  if (!palette %in% c("gene", "protein", "gene_protein", "cell", "compound")) stop('palette should be one of "gene", "protein", "gene_protein", "cell", or "compound".')

  pal <- ddh_pal_c(palette = palette, reverse = reverse)

  ggplot2::scale_color_gradientn(colours = pal(256), ...)
}



#' SCALE COLOR DDH D
#'
#' Color scale constructor for discrete DDH colors
#'
#' @param palette Character name of palette in ddh_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param shuffle Boolean indicating whether the colors should be shuffled
#' @param seed Boolean setting a seed to ensure consistent randomization
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE
#'            or FALSE
#'
#' @examples
#' library(ggplot2)
#' ggplot(mpg, aes(displ, cty, color = manufacturer)) +
#'   geom_point(size = 4) + scale_color_ddh_d()
#' ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Species)) +
#'   geom_point(size = 4) + scale_color_ddh_d("protein", shuffle = TRUE, seed = 123L)
#'
#' @return A color scale for discrete color palettes for various categories
#'
#' @author Matthew Hirschey & Pol Castellano
#'
#' @export
scale_color_ddh_d <- function(palette = "gene", reverse = FALSE, shuffle = FALSE, seed = NULL, ...) {

  palette <- stringr::str_to_lower(palette)

  if (!palette %in% c("gene", "protein", "gene_protein", "cell", "compound")) stop('palette should be one of "gene", "protein", "gene_protein", "cell", or "compound".')
  if (!is.logical(reverse)) stop('reverse should be logical.')
  if (!is.logical(shuffle)) stop('shuffle should be logical.')

  pal <- ddh_pal_d(palette = palette, reverse = reverse, shuffle = shuffle, seed = seed)

  ggplot2::discrete_scale("colour", paste0("ddh_", palette), palette = pal, ...)
}



#' SCALE FILL DDH C
#'
#' FillS scale constructor for continuous DDH color palettes
#'
#' @param palette Character name of palette in ddh_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_fill_gradientn(), used respectively when discrete is TRUE
#'            or FALSE
#'
#' @examples
#' library(ggplot2)
#' ggplot(iris, aes(Sepal.Width, Sepal.Length, fill = Sepal.Width)) +
#'   geom_point(size = 4, shape = 21) + scale_fill_ddh_c()
#' ggplot(iris, aes(Sepal.Width, Sepal.Length, fill = Sepal.Width)) +
#'   geom_point(size = 4, shape = 21) +
#'   scale_fill_ddh_c("Gene_Protein", reverse = TRUE) + theme_ddh()
#'
#' @author Matthew Hirschey & Pol Castellano
#'
#' @export
scale_fill_ddh_c <- function(palette = "gene", reverse = FALSE, ...) {

  palette <- stringr::str_to_lower(palette)

  if (!palette %in% c("gene", "protein", "gene_protein", "cell", "compound")) stop('palette should be one of "gene", "protein", "gene_protein", "cell", or "compound".')

  pal <- ddh_pal_c(palette = palette, reverse = reverse)

  ggplot2::scale_fill_gradientn(colours = pal(256), ...)
}



#' SCALE FILL DDH D
#'
#' FillS scale constructor for discrete DDH color palettes
#'
#' @param palette Character name of palette in ddh_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param shuffle Boolean indicating whether the colors should be shuffled
#' @param seed Boolean setting a seed to ensure consistent randomization
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_fill_gradientn(), used respectively when discrete is TRUE
#'            or FALSE
#'
#' @examples
#' library(ggplot2)
#' ggplot(mpg, aes(class, fill = class)) + geom_bar() + scale_fill_ddh_d()
#' ggplot(mpg, aes(class, fill = class)) + geom_bar() +
#'   scale_fill_ddh_d("protein", reverse = TRUE)
#' ggplot(mpg, aes(class, fill = class)) + geom_bar() +
#'   scale_fill_ddh_d("cell", shuffle = TRUE, seed = 1L)
#'
#' @author Matthew Hirschey & Pol Castellano
#'
#' @export
scale_fill_ddh_d <- function(palette = "gene", reverse = FALSE, shuffle = FALSE, seed = NULL, ...) {

  palette <- stringr::str_to_lower(palette)

  if (!palette %in% c("gene", "protein", "gene_protein", "cell", "compound")) stop('palette should be one of "gene", "protein", "gene_protein", "cell", or "compound".')
  if (!is.logical(reverse)) stop('reverse should be logical.')
  if (!is.logical(shuffle)) stop('shuffle should be logical.')

  pal <- ddh_pal_d(palette = palette, reverse = reverse, shuffle = shuffle, seed = seed)

  ggplot2::discrete_scale("fill", paste0("ddh_", palette), palette = pal, ...)
}

