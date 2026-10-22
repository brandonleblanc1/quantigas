#' Re-export [magrittr::%>%]
#'
#' @importFrom magrittr %>%
#'
#' @export
magrittr::`%>%`

#' [ggplot2::ggplot] Publication Theme
#'
#' Add a sleek publication-style theme to a [ggplot2::ggplot] object.
#'
#' Copied and modified from https://github.com/koundy/ggplot_theme_Publication.
#'
#' @param base_size (`numeric(1)`)
#'   Font size used in the plot. Defaults to 11.
#' @param base_family (`character(1)`)
#'   Font family to use in the plot. Defaults to "sans".
#'
#' @export
theme_publication <- function(
  base_size = defaults$FONTSIZE_DEFAULT,
  base_family = defaults$FONTFAMILY_DEFAULT
) {
  (ggthemes::theme_foundation(
    base_size = base_size,
    base_family = base_family
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(1.2),
        hjust = 0.5,
        margin = ggplot2::margin(0, 0, 20, 0)
      ),
      text = ggplot2::element_text(),
      panel.background = ggplot2::element_rect(color = "white", fill = "white"),
      plot.background = ggplot2::element_rect(color = NA),
      panel.border = ggplot2::element_rect(color = NA),
      axis.title = ggplot2::element_text(
        face = "plain",
        size = ggplot2::rel(1)
      ),
      axis.title.y = ggplot2::element_text(angle = 90, vjust = 2),
      axis.title.x = ggplot2::element_text(vjust = -0.2),
      axis.text = ggplot2::element_text(),
      axis.line.x = ggplot2::element_line(color = "black"),
      axis.line.y = ggplot2::element_line(color = "black"),
      axis.ticks = ggplot2::element_line(),
      panel.grid.major = ggplot2::element_line(color = "#f0f0f0"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.key = ggplot2::element_rect(color = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.key.size = grid::unit(0.5, "cm"),
      legend.title = ggplot2::element_text(face = "italic"),
      plot.margin = grid::unit(c(10, 5, 5, 5), "mm"),
      strip.background = ggplot2::element_rect(
        color = "#f0f0f0",
        fill = "#f0f0f0"
      ),
      strip.text = ggplot2::element_text(face = "bold")
    ))
}
