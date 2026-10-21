#' [R6::R6Class] for Model Building
#'
#' Derived class for modelling GSMG measurement data. Specifically deals with
#' building an error-in-variables model using polynomials for the analysis
#' function. It inherits from the [ModelBuilder] class and thus has access to
#' its public and private attributes and methods.
#'
#' @export
PolyModel <- R6::R6Class(
  classname = "PolyModel",
  lock_class = TRUE,
  inherit = ModelBuilder,
  public = list(
    #' @description
    #' Print Object
    #'
    #' @return The calling object invisibly.
    print = function() {
      super$print()

      cat(
        "<PolyModel>\n",
        "  Private Attributes:\n",
        "    ratio mc samples: matrix [",
        paste0(dim(private$r_mc), collapse = " x "),
        "]\n",
        "    amount mc samples: matrix [",
        paste0(dim(private$c_mc), collapse = " x "),
        "]\n",
        "    solution mc samples: matrix [",
        paste0(dim(private$solution_mc), collapse = " x "),
        "]\n",
        "    mc residuals summary: tbl_df [",
        paste0(dim(private$mc_summ), collapse = " x "),
        "]\n",
        "    coefficients data: tbl_df [",
        paste0(dim(private$beta_cov), collapse = " x "),
        "]\n",
        "    ratio residuals plot: ",
        ifelse(is.null(private$r_resid_plot), "NULL\n", "ggplot\n"),
        "    amount residuals plot: ",
        ifelse(is.null(private$c_resid_plot), "NULL\n", "ggplot\n"),
        "    ratio mc residuals plot: ",
        ifelse(is.null(private$r_mc_plot), "NULL\n", "ggplot\n"),
        "    amount mc residuals plot: ",
        ifelse(is.null(private$c_mc_plot), "NULL\n", "ggplot\n"),
        "    analysis function plot: ",
        ifelse(is.null(private$model_plot), "NULL\n", "ggplot\n"),
        sep = ""
      )

      invisible(self)
    },

    #' @description
    #' Plot Selected Model
    #'
    #' @details
    #' This function generates a plot of the selected model's analysis function
    #' that additionally displays information relating to the fit of the
    #' selected model.
    #'
    #' @param model_num (`numeric(1)`)
    #'   Number identifying the candidate polynomial model to plot. Defaults to
    #'   1.
    #' @param units (`character(1)`)
    #'   Units to use for the y-axis label. Defaults to "mu*mol ~~ mol^-1".
    #' @param fontfamily (`character(1)`)
    #'   Font family to use for plot text. Defaults to "sans".
    #'
    #' @return The calling object invisibly.
    plot_model = function(
      model_num = defaults$MODEL_NUM,
      units = defaults$UNITS,
      fontfamily = defaults$FONTFAMILY
    ) {
      model_info <- private$models %>%
        dplyr::filter(.data$model_num == .env$model_num)

      beta <- model_info$beta[[1]]

      # Building model formula to display on plot.
      model_form <- c(
        "c ==",
        purrr::imap_chr(signif(beta, 5), \(val, i) {
          if (i == 1) {
            as.character(val)
          } else {
            coef <- if (val > 0) "+ {val}" else "- {-val}"
            stringr::str_glue("{stringr::str_glue(coef)} ~~ r^{i - 1}")
          }
        })
      ) %>%
        stringr::str_flatten(collapse = " ")

      # Generating text elements to display on plot.
      plot_text <- if (inherits(self, "PolyModel_Replicates")) {
        list(
          grid::textGrob(
            str2expression(paste("R^2 ==", signif(model_info$adj_r_sq))),
            x = 0.15,
            y = 0.9,
            just = "left"
          ),
          grid::textGrob(
            str2expression(paste("tau ==", signif(model_info$tau))),
            x = 0.15,
            y = 0.85,
            just = "left"
          ),
          grid::textGrob(
            str2expression(paste("nu ==", signif(model_info$nu, digits = 3))),
            x = 0.15,
            y = 0.8,
            just = "left"
          )
        )
      }

      # Generating points defined by the analysis function.
      r_vals <- seq(
        min(private$cal_data$r),
        max(private$cal_data$r),
        length = 1000
      )
      c_vals <- vapply(
        r_vals,
        \(r) {
          eval(str2expression(model_info$funcs))(r, beta)
        },
        numeric(1)
      )

      private$model_plot <- private$cal_data %>%
        ggplot2::ggplot(ggplot2::aes(x = r, y = c)) +
        ggplot2::geom_line(
          data = tibble::tibble(r = r_vals, c = c_vals),
          alpha = 0.5,
          linetype = 2
        ) +
        ggplot2::geom_pointrange(
          ggplot2::aes(ymin = c - uc, ymax = c + uc),
          size = 0.25
        ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(xmin = r - ur, xmax = r + ur),
          width = 0
        ) +
        ggplot2::annotation_custom(rlang::inject(grid::grobTree(
          grid::textGrob(
            str2expression(model_form),
            x = 0.15,
            y = 0.95,
            just = "left"
          ),
          !!!plot_text,
          gp = grid::gpar(
            col = "#494F55",
            fontsize = 11,
            fontfamily = fontfamily
          )
        ))) +
        ggplot2::labs(
          x = "Ratio",
          y = str2expression(paste0("Amount ~~ Fraction ~~ (", units, ")"))
        ) +
        ggplot2::theme_bw()

      invisible(self)
    },

    #' @description
    #' Get MC Samples
    #'
    #' @details
    #' This function returns a list of all private attributes containing MC
    #' samples. If all attributes are not needed at the same time, it can be
    #' convenient to select a specific one from the function call, i.e.
    #' `get_mc()$r_mc`.
    #'
    #' @return A `list` of attributes.
    get_mc = function() {
      list(
        r_mc = private$r_mc,
        c_mc = private$c_mc,
        solution_mc = private$solution_mc
      )
    },

    #' @description
    #' Get Data Frames
    #'
    #' @details
    #' This function returns a list of all private attributes that are data
    #' frames. If all attributes are not needed at the same time, it can be
    #' convenient to select a specific one from the function call, i.e.
    #' `get_data()$cal_data`.
    #'
    #' @return A `list` of attributes.
    get_data = function() {
      rlang::list2(
        !!!super$get_data(),
        mc_summ = private$mc_summ,
        beta_cov = private$beta_cov
      )
    },

    #' @description
    #' Get Plots
    #'
    #' @details
    #' This function returns a list of all private attributes that are plots.
    #' If all attributes are not needed at the same time, it can be convenient
    #' to select a specific one from the function call, i.e.
    #' `get_plots()$r_resid_plot`.
    #'
    #' @return A `list` of attributes.
    get_plots = function() {
      list(
        r_resid_plot = private$r_resid_plot,
        c_resid_plot = private$c_resid_plot,
        r_mc_plot = private$r_mc_plot,
        c_mc_plot = private$c_mc_plot,
        model_plot = private$model_plot
      )
    }
  ),
  private = list(
    # `matrix`: Monte Carlo samples for the replicate ratios.
    r_mc = NULL,

    # `matrix`: Monte Carlo samples for the amount fractions.
    c_mc = NULL,

    # `matrix`: Monte Carlo samples for the optimized model parameters.
    solution_mc = NULL,

    # [tibble::tibble]: Summary of residuals from the Monte Carlo simulation
    # for both ratios and amount fractions.
    mc_summ = NULL,

    # [tibble::tibble]: Covariance of coefficients for chosen model.
    beta_cov = NULL,

    # [ggplot2::ggplot]: Plot of the standardized ratio residuals for each
    # candidate model.
    r_resid_plot = NULL,

    # [ggplot2::ggplot]: Plot of the standardized amount fraction residuals
    # for each candidate model.
    c_resid_plot = NULL,

    # [ggplot2::ggplot]: Plot of the ratio residuals from the MC simulation.
    r_mc_plot = NULL,

    # [ggplot2::ggplot]: Plot of the amount fraction residuals from the MC
    # simulation.
    c_mc_plot = NULL,

    # [ggplot2::ggplot]: Plot of the selected model's analysis function.
    model_plot = NULL,

    # This function builds the formulas and analysis functions for candidate
    # polynomial models up to the 3rd degree.
    build_models = function() {
      degrees <- seq(0, 3)

      # Building formulas to use in `lm()`.
      forms <- degrees %>%
        purrr::map_chr(\(num) {
          if (num == 0) "c ~ 1" else paste0("I(r^", num, ")")
        }) %>%
        purrr::accumulate(
          \(form1, form2) paste(form1, form2, sep = " + "),
          .simplify = TRUE
        )

      # Building analysis functions.
      funcs <- degrees %>%
        purrr::accumulate(c) %>%
        purrr::map_chr(\(set) {
          paste0(
            "\\(rho, beta) sum(beta * rho^c(",
            paste0(set, collapse = ", "),
            "))"
          )
        })

      # Removing 0th degree elements.
      tibble::tibble(
        model_num = degrees[-1],
        forms = forms[-1],
        funcs = funcs[-1]
      )
    },

    # This function generates a plot of the standardized residuals for each
    # candidate polynomial model for either the ratios or the amount fractions.
    plot_models_resid = function(plot_var = c("r", "c")) {
      plot_var <- match.arg(plot_var)
      fitted_var <- switch(plot_var, r = "rho", c = "xi")
      zscore <- paste0(plot_var, "_zscore")

      private$models %>%
        dplyr::select(model_num, tidyselect::all_of(c(fitted_var, zscore))) %>%
        {
          if (inherits(self, "PolyModel_Replicates") && (plot_var == "r")) {
            dplyr::mutate(
              .,
              rho = purrr::map(rho, \(ratios) rep(ratios, private$counts$n))
            )
          } else {
            .
          }
        } %>%
        tidyr::unnest_longer(!model_num) %>%
        ggplot2::ggplot(ggplot2::aes(
          x = .data[[fitted_var]],
          y = .data[[zscore]]
        )) +
        ggplot2::geom_hline(
          ggplot2::aes(yintercept = 2),
          linetype = 2,
          color = "red"
        ) +
        ggplot2::geom_hline(
          ggplot2::aes(yintercept = -2),
          linetype = 2,
          color = "red"
        ) +
        ggplot2::geom_smooth(
          method = "loess",
          formula = y ~ x,
          se = FALSE,
          linetype = 3,
          color = "black",
          alpha = 0.5
        ) +
        ggplot2::geom_point(size = 1.5, color = "black") +
        ggplot2::facet_wrap(
          ggplot2::vars(model_num),
          ncol = 2,
          scales = "free",
          labeller = "label_both"
        ) +
        ggplot2::labs(
          x = stringr::str_to_title(fitted_var),
          y = "Standardized Residual"
        ) +
        ggplot2::theme_bw()
    },

    # This function generates a plot of the selected model's MC sample
    # residuals for either the ratios or the amount fractions.
    plot_mc_resid = function(plot_var = c("r", "c"), mc_data) {
      plot_var <- match.arg(plot_var)

      tlabel <- switch(
        plot_var,
        r = "Ratio Residuals",
        c = "Amount Fraction Residuals"
      )

      resid <- paste0(plot_var, "_resid")
      uresid <- paste0(plot_var, "_uresid")

      # Getting unnested residuals for violin plots.
      unnested <- mc_data %>%
        tibble::as_tibble(
          .name_repair = "unique_quiet",
          rownames = "sample_id"
        ) %>%
        tidyr::pivot_longer(tidyselect::starts_with("..."))

      ggplot2::ggplot(
        mapping = ggplot2::aes(x = forcats::as_factor(sample_id))
      ) +
        ggplot2::geom_hline(
          ggplot2::aes(yintercept = 0),
          linetype = 2,
          color = "red",
          alpha = 0.5
        ) +
        ggplot2::geom_violin(
          ggplot2::aes(y = value),
          data = unnested,
          quantile.linetype = 3,
          quantile.color = "grey",
          alpha = 0.8
        ) +
        ggplot2::geom_pointrange(
          ggplot2::aes(
            y = .data[[resid]],
            ymin = .data[[resid]] - .data[[uresid]],
            ymax = .data[[resid]] + .data[[uresid]]
          ),
          data = private$mc_summ
        ) +
        ggplot2::labs(title = tlabel, x = NULL, y = "Residual") +
        ggplot2::theme_bw()
    },

    # This function calculates the covariances and correlations of the model
    # coefficients.
    calc_beta_cov = function(beta, beta_mc) {
      n_coef <- length(beta)
      row_names <- seq_along(beta) %>%
        purrr::map_chr(\(i) stringr::str_c("b", i))
      col_names <- c("value", "uncert", rep("cov", n_coef), rep("cor", n_coef))

      robustbase::covMcd(t(beta_mc))$cov %>%
        rbind(beta, sqrt(diag(.)), ., stats::cov2cor(.)) %>%
        t() %>%
        `rownames<-`(row_names) %>%
        `colnames<-`(col_names) %>%
        tibble::as_tibble(
          .name_repair = "unique_quiet",
          rownames = "coefficient"
        ) %>%
        dplyr::relocate(coefficient)
    }
  )
)
