#' Function to generate color palettes
#'
#' @param hex Hex code for main color as string
#' @param discrete Reorder for categorical data to place main color first
#' @return A set of three colors
#'
#' @examples
#' generate_colors("#2EC09C")
#' generate_colors("#BE34EF")
generate_colors <- function(hex, discrete = FALSE) {
  color_set <- c(
    colorspace::desaturate(colorspace::lighten(hex, .6), .2),
    hex,
    colorspace::darken(hex, .6, space = "HLS")
  )

  if (discrete == TRUE) color_set <- color_set[c(2, 1, 3)]

  return(color_set)
}


#' Function to generate shuffled colors based on `colorRampPalette`
#'
#' @param colors Colors to interpolate; must be a valid argument to col2rgb().
#' @return A function with argument a vector of values between 0 and 1 that are
#'         mapped to a numeric matrix of RGB color values with one row per color
#'         and 3 or 4 columns
#'
#' @examples
#' colorRampPalette_shuffle(colors = c("#DEBBF4", "#BE34EF", "#AC7B84"))
#' colorRampPalette_shuffle(colors = c("#DEBBF4", "#BE34EF", "#AC7B84"), seed = 123L)
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

#' Function to extract DDH colors as hex codes
#'
#' @param ... Character names of colors
#'
#' @examples
#' ddh_colors()
#' ddh_colors("gene")
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


#' Return function to interpolate a continuous DDH color palette
#'
#' @param palette Character name of palette in ddh_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments to pass to [grDevices::colorRampPalette()]
#'
#' @examples
#' ddh_pal_c()(10)
#' ddh_pal_c("protein")(3)
#' ddh_pal_c("cOMPOUND", reverse = TRUE)(5)
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


#' Return function to interpolate a discrete DDH color palette
#'
#' @param palette Character name of palette in ddh_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#'
#' @examples
#' ddh_pal_d()(10)
#' ddh_pal_d("protein")(3)
#' ddh_pal_d("cOMPOUND", reverse = TRUE)(5)
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


#' Color scale constructor for continuous DDH color palettes
#'
#' @param palette Character name of palette in ddh_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE
#'            or FALSE
#' @examples
#' library(ggplot2)
#' ggplot(mpg, aes(hwy, cty, color = displ)) + geom_point(size = 4) +
#' scale_color_ddh_c()
#' ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Sepal.Width)) +
#'   geom_point(size = 4) + scale_color_ddh_c("cOMPOUND", reverse = TRUE)
#' @export
scale_color_ddh_c <- function(palette = "gene", reverse = FALSE, ...) {

  palette <- stringr::str_to_lower(palette)

  if (!palette %in% c("gene", "protein", "gene_protein", "cell", "compound")) stop('palette should be one of "gene", "protein", "gene_protein", "cell", or "compound".')

  pal <- ddh_pal_c(palette = palette, reverse = reverse)

  ggplot2::scale_color_gradientn(colours = pal(256), ...)
}



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
#' @export
scale_color_ddh_d <- function(palette = "gene", reverse = FALSE, shuffle = FALSE, seed = NULL, ...) {

  palette <- stringr::str_to_lower(palette)

  if (!palette %in% c("gene", "protein", "gene_protein", "cell", "compound")) stop('palette should be one of "gene", "protein", "gene_protein", "cell", or "compound".')
  if (!is.logical(reverse)) stop('reverse should be logical.')
  if (!is.logical(shuffle)) stop('shuffle should be logical.')

  pal <- ddh_pal_d(palette = palette, reverse = reverse, shuffle = shuffle, seed = seed)

  ggplot2::discrete_scale("colour", paste0("ddh_", palette), palette = pal, ...)
}



#' Fill scale constructor for continuous DDH color palettes
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
#' @export
scale_fill_ddh_c <- function(palette = "gene", reverse = FALSE, ...) {

  palette <- stringr::str_to_lower(palette)

  if (!palette %in% c("gene", "protein", "gene_protein", "cell", "compound")) stop('palette should be one of "gene", "protein", "gene_protein", "cell", or "compound".')

  pal <- ddh_pal_c(palette = palette, reverse = reverse)

  ggplot2::scale_fill_gradientn(colours = pal(256), ...)
}



#' Fill scale constructor for discrete DDH color palettes
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
#' @export
scale_fill_ddh_d <- function(palette = "gene", reverse = FALSE, shuffle = FALSE, seed = NULL, ...) {

  palette <- stringr::str_to_lower(palette)

  if (!palette %in% c("gene", "protein", "gene_protein", "cell", "compound")) stop('palette should be one of "gene", "protein", "gene_protein", "cell", or "compound".')
  if (!is.logical(reverse)) stop('reverse should be logical.')
  if (!is.logical(shuffle)) stop('shuffle should be logical.')

  pal <- ddh_pal_d(palette = palette, reverse = reverse, shuffle = shuffle, seed = seed)

  ggplot2::discrete_scale("fill", paste0("ddh_", palette), palette = pal, ...)
}

