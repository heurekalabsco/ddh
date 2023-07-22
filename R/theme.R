#' THEME_DDH
#'
#' A ggplot theme which allow custom yet consistent styling of plots in the
#' Data-driven hypothesis web app.
#'
#' @param style (string) Overall color style of text labels.
#' Either "default" or "light".
#' @param base_size (integer) Base point size
#' @param grid (string) Grid lines. Options include  "none" or
#' any combination of "X", "Y", "x" and "y".
#' @param axistitle (string) Axis titles. Options include "none" or
#' any combination of "X", "Y", "x" and "y".
#' @param axistext (string) Axis text labels for values or groups.
#' Options include "none" or any combination of "X", "Y", "x" and "y".
#' @param axis_num (string) Should axis text be formatted as monospaced? Set
#' to  x and y, respectively, in case numeric values are displayed. Options
#' include "none" or any combination of "X", "Y", "x" and "y".
#' @param legend_num (logical) Should legend text be formatted as monospaced?
#' Defaults to FALSE (no monospace font). Set to TRUE in case of numeric values.
#' @param panelborder (logical) Should a panel border be drawn?
#' Defaults to FALSE (no border). If TRUE it also adds tick marks to both axes.
#' @param margin (numeric) Should a margin of x be added to the plot?
#' Defaults to 0 (no margin by default).
#' @param ... Other arguments passed to ggplot methods.
#'
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(diamonds, aes(cut)) + geom_bar() + theme_ddh()
#' ggplot(diamonds, aes(cut)) + geom_bar() +
#'   theme_ddh(style = "light", axisline = "none")
#' ggplot(diamonds, aes(cut)) + geom_bar() +
#'   theme_ddh(grid = "y", axisline = "none", panelborder = TRUE,
#'             panelbackground = TRUE, axis_num = "y", margin = 20)
#' }
#'
#' @author Matthew Hirschey & Pol Castellano
#'
#' @export
theme_ddh <- function(style = "default", base_size = 15,
                      grid = "xy", axisline = "xy",
                      axistitle = "xy", axistext = "xy",
                      axis_num = "none", legend_num = FALSE,
                      panelborder = FALSE, panelbackground = FALSE,
                      margin = 2, ...) {
  if(!style %in% c("default", "light")) stop('style must be either "default" or "light"')
  if(!is.character(grid)) stop('grid must be a character: "none" or any combination of "X", "Y", "x" and "y"')
  if(!is.character(axisline)) stop('axisline must be a character: "none" or any combination of "X", "Y", "x" and "y"')
  if(!is.character(axistitle)) stop('axistitle must be a character: "none" or any combination of "X", "Y", "x" and "y"')
  if(!is.character(axistext)) stop('axistext must be a character: "none" or any combination of "X", "Y", "x" and "y"')
  if(!is.character(axis_num)) stop('axis_num must be a character: "none" or any combination of "X", "Y", "x" and "y"')
  if(!is.logical(legend_num)) stop('legend_num must be a logical variable')
  if(!is.logical(panelborder)) stop('panelborder must be a logical variable')
  if(!is.logical(panelbackground)) stop('panelbackground must be a logical variable')
  if(!is.numeric(margin)) stop('margin must be a numeric value')

  fontfamily_sans <- "Nunito Sans"
  fontfamily_slab <- "Roboto Slab"
  fontfamily_mono <- "Roboto Mono"

  if (style == "default") { base_col <- "black"; light_col <- "grey20" }
  if (style == "light") { base_col <- "grey25"; light_col <- "grey45" }

  out <-
    ggplot2::theme_minimal(base_size = base_size, base_family = fontfamily_slab) +
    ggplot2::theme(
      text = ggplot2::element_text(
        color = base_col
      ),
      line = ggplot2::element_line(
        color = light_col
      ),
      rect = ggplot2::element_rect(
        color = light_col,
        fill = "transparent"
      ),
      plot.title = ggtext::element_textbox_simple(
        family = fontfamily_sans,
        size = base_size * 1.7,
        face = "bold",
        lineheight = .8,
        box.color = NA,
        margin = ggplot2::margin(t = 0, b = base_size * .67)
      ),
      plot.subtitle = ggtext::element_textbox_simple(
        family = fontfamily_sans,
        size = base_size,
        lineheight = 1.2,
        color = light_col,
        margin = ggplot2::margin(t = 0, b = base_size * 1.5)
      ),
      plot.caption = ggtext::element_textbox_simple(
        family = fontfamily_sans,
        size = base_size / 2,
        lineheight = 1.2,
        color = light_col,
        margin = ggplot2::margin(t = base_size * 1.5, b = 0)
      ),
      axis.title.x = ggplot2::element_text(
        family = fontfamily_sans,
        margin = ggplot2::margin(t = base_size / 3, r = 3, b = 3, l = 3)
      ),
      axis.title.y = ggplot2::element_text(
        family = fontfamily_sans,
        margin = ggplot2::margin(t = 3, r = 3, b = base_size / 3, l = 3)
      ),
      axis.text.x = ggplot2::element_text(
        color = light_col,
        margin = ggplot2::margin(t = base_size / 4, r = 1, b = 1, l = 1)
      ),
      axis.text.y = ggplot2::element_text(
        color = light_col,
        margin = ggplot2::margin(t = 1, r = base_size / 4, b = 1, l = 1)
      ),
      axis.ticks.length = grid::unit(.33, "lines"),
      strip.text = ggplot2::element_text(
        family = fontfamily_sans,
        face = "bold"
      ),
      strip.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = "grey87",
        size = .4
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(
        color = "transparent",
        fill = "transparent"
      ),
      plot.background = ggplot2::element_rect(
        color = "white",
        fill = "white"
      ),
      plot.margin = ggplot2::margin(
        t = margin, r = margin, l = margin, b = margin
      ),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      legend.position = "right",
      legend.title = ggplot2::element_text(
        family = fontfamily_sans,
        color = light_col,
        size = base_size * .75,
        margin = ggplot2::margin(b = 10)
      ),
      legend.text = ggplot2::element_text(
        family = fontfamily_slab,
        color = light_col,
        size = base_size * .75
      )
    )

  ## add panel border if selected
  if (panelborder == TRUE) {
    out <- out +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(
          color = light_col,
          fill = "transparent",
          size = .65
        ),
        axis.ticks = ggplot2::element_line(
          color = light_col,
          size = .4
        )
      )
  }

  ## add panel background if selected
  if (panelbackground == TRUE) {
    out <- out +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(
          color = "transparent",
          fill = "grey95",
          size = .5
        ),
        panel.grid.major = ggplot2::element_line(
          color = "grey82",
          size = .4
        )
      )
  }

  ## add panel grid if selected
  if (grid != "none") {
    if (!stringr::str_detect(grid, "X|x")) {
      out <- out +
        ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
    }
    if (!stringr::str_detect(grid, "Y|y")) {
      out <- out +
        ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
    }
  } else {
    out <- out + ggplot2::theme(panel.grid.major = ggplot2::element_blank())
  }

  ## add axis line(s) if selected
  if (axisline != "none") {
    if (stringr::str_detect(axisline, "X|x")) {
      out <- out +
        ggplot2::theme(
          axis.line.x = ggplot2::element_line(color = light_col, size = .65),
          axis.ticks.x = ggplot2::element_line(color = light_col, size = .4),
          axis.text.x = ggplot2::element_text(
            margin = ggplot2::margin(t = base_size / 2, r = 0, b = 0, l = 0)
          ))
    }
    if (stringr::str_detect(axisline, "Y|y")) {
      out <- out +
        ggplot2::theme(
          axis.line.y = ggplot2::element_line(color = light_col, size = .65),
          axis.ticks.y = ggplot2::element_line(color = light_col, size = .4),
          axis.text.y = ggplot2::element_text(
            margin = ggplot2::margin(t = 0, r = base_size / 2, b = 0, l = 0)
          ))
    }
  }

  ## remove axis titles if selected
  if (axistitle != "none") {
    if (!stringr::str_detect(axistitle, "X|x")) {
      out <- out +
        ggplot2::theme(axis.title.x = ggplot2::element_blank())
    }
    if (!stringr::str_detect(axistitle, "Y|y")) {
      out <- out +
        ggplot2::theme(axis.title.y = ggplot2::element_blank())
    }
  } else {
    out <- out +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank())
  }

  ## turn axis text into monospaced if selected
  if (axis_num != "none") {
    if (stringr::str_detect(axis_num, "X|x")) {
      out <- out +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
          family = fontfamily_mono,
          color = light_col,
          margin = ggplot2::margin(t = base_size / 4, r = 1, b = 1, l = 1)
        ))
    }
    if (stringr::str_detect(axis_num, "Y|y")) {
      out <- out +
        ggplot2::theme(axis.text.y = ggplot2::element_text(
          family = fontfamily_mono,
          color = light_col,
          margin = ggplot2::margin(t = base_size / 4, r = 1, b = 1, l = 1)
        ))
    }
  }

  ## remove axis text if selected
  if (axistext != "none") {
    if (!stringr::str_detect(axistext, "X|x")) {
      out <- out +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank())
    }
    if (!stringr::str_detect(axistext, "Y|y")) {
      out <- out +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())
    }
  } else {
    out <- out +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())
  }

  ## turn legend text into monospaced if selected
  if (legend_num == TRUE) {
    out <- out +
      ggplot2::theme(legend.text = ggplot2::element_text(
        family = fontfamily_mono,
        color = light_col,
        size = base_size * .75,
        hjust = 1 ## right-adjusted to allow for correct alignment
      ))
  }

  return(out)
}
