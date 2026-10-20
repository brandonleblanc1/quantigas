#' @title
#' ggplot2 Publication Theme
#'
#' @description
#' Add a sleek publication-style theme to a ggplot2 object.
#'
#' @details
#' Copied and modified from https://github.com/koundy/ggplot_theme_Publication.
#'
#' @param base_size (`numeric(1)`)
#'   Font size used in the plot. Defaults to 11.
#' @param base_family (`character(1)`)
#'   Font family to use in the plot. Defaults to "sans".
#'
#' @examples
#' plot_to_customize + theme_Publication()
theme_publication <- function(base_size = 11, base_family = "sans") {
  (ggthemes::theme_foundation(
    base_size = base_size,
    base_family = base_family
  ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = rel(1.2),
        hjust = 0.5,
        margin = margin(0, 0, 20, 0)
      ),
      text = element_text(),
      panel.background = element_rect(color = "white", fill = "white"),
      plot.background = element_rect(color = NA),
      panel.border = element_rect(color = NA),
      axis.title = element_text(face = "plain", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(color = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(color = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.key.size = unit(0.5, "cm"),
      legend.title = element_text(face = "italic"),
      plot.margin = unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(color = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    ))
}
