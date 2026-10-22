#' [R6::R6Class] for Sample Data
#'
#' Derived class for processing data from COGAS. Specifically deals with data
#' for unknown samples and calculating/visualizing results from the calibration
#' model. Inherits from the [COGASData] class and thus has access to its public
#' and private attributes and methods.
#'
#' @export
SampData <- R6::R6Class(
  classname = "SampData",
  lock_class = TRUE,
  inherit = COGASData,
  public = list(
    #' @description
    #' Create a new `SampData` object.
    #'
    #' @param params (`list`)
    #'   Parameters to initialize a [COGASData] object. The private attributes
    #'   can either be given as normal, or a [CalData] object can be specified
    #'   to pull information from.
    #'
    #' @return A new `SampData` object.
    initialize = function(params) {
      # Validating inputs.
      stopifnot(length(params) > 0)

      # Allowing two different ways to initialize a `SampData` object: 1) pass
      # a [CalData] object and copy its parameters, or 2) treat the
      # initialization as for a normal [COGASData] object.
      if ("cal_obj" %in% names(params)) {
        stopifnot(inherits(params$cal_obj, "CalData"))

        cal_info <- cal_obj$get_info()

        if (inherits(self, "SampData_Embedded")) {
          super$initialize(rlang::list2(
            !!!cal_info,
            filepath = NULL,
            !!!params
          ))
        } else {
          super$initialize(rlang::list2(
            !!!cal_info,
            !!!params[names(params) != "cal_obj"]
          ))
        }
      } else {
        super$initialize(params)
      }
    },

    #' @description
    #' Print useful information about the state of the object.
    #'
    #' @return The calling object invisibly.
    print = function() {
      super$print()

      cat(
        "<SampData>\n",
        "  Private Attributes:\n",
        "    amount mc samples: matrix [",
        paste0(dim(private$amount_mc), collapse = " x "),
        "]\n",
        "    calculated amounts: tbl_df [",
        paste0(dim(private$amounts), collapse = " x "),
        "]\n",
        "    amount violin plot: ",
        ifelse(is.null(private$amount_violin_plot), "NULL\n", "ggplot\n"),
        sep = ""
      )

      invisible(self)
    },

    #' @description
    #' Calculate amount fractions and uncertainties for unknown samples based
    #' on either a calibration model for multi-point calibration or through a
    #' one-to-one comparison with the control sample for single-point
    #' calibration. The final amount fractions are estimated from the ratio MC
    #' samples that were used in the ratio uncertainty calculation.
    #'
    #' @param model_obj ([ModelBuilder] | `numeric` | [tibble::tibble])
    #'   Calibration model object, or a vector/list/data frame containing the
    #'   following named parameters for the control sample:
    #'
    #'   * amount: certified amount fraction
    #'   * uamount: uncertainty associated with the amount fraction
    #'
    #' @param summ_cols (`character`)
    #'   Variables to group by when performing calculation. Defaults to
    #'   "sample_id".
    #'
    #' @return The calling object invisibly.
    calc_amounts = function(model_obj, summ_cols = defaults$SUMM_COLS) {
      if (inherits(model_obj, "ModelBuilder")) {
        control_mc <- model_obj$predict()
        private$amount_mc <- model_obj$predict(private$ratio_mc)
      } else {
        stopifnot(c("amount", "uamount") %in% names(model_obj))

        # Generating MC samples for the control.
        control_mc <- stats::rnorm(
          ncol(private$ratio_mc),
          mean = model_obj[["amount"]],
          sd = model_obj[["uamount"]]
        )
        private$amount_mc <- sweep(private$ratio_mc, 2, control_mc, FUN = "*")
      }

      # Summarizing results.
      private$amounts <- private$amount_mc %>%
        tibble::as_tibble(
          .name_repair = "unique_quiet",
          rownames = "sample_set"
        ) %>%
        dplyr::left_join(private$cleaned_data, ., by = "sample_set") %>%
        tidyr::pivot_longer(tidyselect::starts_with("...")) %>%
        dplyr::summarise(
          amount = robustbase::huberM(value)$mu,
          uamount = robustbase::Qn(value),
          .by = tidyselect::all_of(summ_cols)
        ) %>%
        tibble::add_row(
          sample_id = as.factor(private$control),
          amount = robustbase::huberM(control_mc)$mu,
          uamount = robustbase::Qn(control_mc),
          .before = 1
        ) %>%
        dplyr::mutate(uamount_rel = 100 * uamount / amount)

      invisible(self)
    },

    #' @description
    #' Generate useful diagnostic plots for the current analysis.
    #'
    #' @param units (`character(1)`)
    #'   Units to use for the y-axis label. Defaults to "mu*mol ~~ mol^-1".
    #'
    #' @return The calling object invisibly.
    plot_data = function(units = defaults$UNITS) {
      if (!inherits(self, "SampData_Embedded")) {
        super$plot_data()
      }

      if (!is.null(private$amounts)) {
        # Getting unnested amounts for violin plots.
        unnested <- private$amount_mc %>%
          tibble::as_tibble(
            .name_repair = "unique_quiet",
            rownames = "sample_set"
          ) %>%
          dplyr::left_join(private$cleaned_data, ., by = "sample_set") %>%
          tidyr::pivot_longer(tidyselect::starts_with("..."))

        # If an additional column was used to nest data during the calculation,
        # the violin plots will get filled accordingly to visualize the
        # different groupings.
        additional_cols <- stringr::str_ends(
          names(private$amounts),
          "(set)|(port)"
        )

        violin_layer <- if (any(additional_cols)) {
          fill_var <- names(private$amounts)[additional_cols][1]
          scale_name <- stringr::str_split(fill_var, "_")[[1]] %>%
            stringr::str_to_title() %>%
            stringr::str_flatten(collapse = " ")

          list(
            ggplot2::geom_violin(
              ggplot2::aes(y = value, fill = .data[[fill_var]]),
              data = unnested,
              quantile.linetype = 3,
              quantile.color = "black",
              scale = "width",
              alpha = 0.6,
              position = "identity"
            ),
            ggplot2::scale_fill_viridis_d(name = scale_name, na.value = "white")
          )
        } else {
          list(ggplot2::geom_violin(
            ggplot2::aes(y = value),
            data = unnested,
            quantile.linetype = 3,
            quantile.color = "grey",
            scale = "width",
            alpha = 0.8
          ))
        }

        # Plotting distribution of amount fractions by sample id, overlaid
        # with their location and uncertainty estimates.
        private$amount_violin_plot <- ggplot2::ggplot(
          mapping = ggplot2::aes(x = sample_id)
        ) +
          violin_layer +
          ggplot2::geom_pointrange(
            ggplot2::aes(
              y = amount,
              ymin = amount - uamount,
              ymax = amount + uamount
            ),
            data = private$amounts,
            size = 0.1
          ) +
          ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2)) +
          ggplot2::labs(
            x = NULL,
            y = str2expression(paste0("Amount ~~ Fraction ~~ (", units, ")"))
          ) +
          ggplot2::theme_bw()
      }

      invisible(self)
    },

    #' @description
    #' Get all private attributes containing MC samples. If all attributes are
    #' not needed at the same time, it can be convenient to select a specific
    #' one from the function call, i.e. `get_mc()$response_mc`.
    #'
    #' @return A `list` of attributes.
    get_mc = function() {
      rlang::list2(!!!super$get_mc(), amount_mc = private$amount_mc)
    },

    #' @description
    #' Get all private attributes that are data frames. If all attributes are
    #' not needed at the same time, it can be convenient to select a specific
    #' one from the function call, i.e. `get_data()$cleaned_data`.
    #'
    #' @return A `list` of attributes.
    get_data = function() {
      rlang::list2(!!!super$get_data(), amounts = private$amounts)
    },

    #' @description
    #' Get all private attributes that are plot objects. If all attributes are
    #' not needed at the same time, it can be convenient to select a specific
    #' one from the function call, i.e. `get_plots()$break_set_plot`.
    #'
    #' @return A `list` of attributes.
    get_plots = function() {
      rlang::list2(
        !!!super$get_plots(),
        amount_violin_plot = private$amount_violin_plot
      )
    }
  ),
  private = list(
    # `matrix`: Monte Carlo samples for calculated amount fractions.
    amount_mc = NULL,

    # [tibble::tibble] Calculated amount fractions and their associated
    # uncertainties based on the calibration model.
    amounts = NULL,

    # [ggplot2::ggplot] Plot of the distributions of the amount fraction MC
    # samples for each sample.
    amount_violin_plot = NULL
  )
)
