#' [R6::R6Class] for Sample Data
#'
#' Derived class for processing data from COGAS. Specifically deals with QA
#' data for unknown samples and comparing a sample's certified value to the
#' current calculation. It inherits from the [SampData] class and thus has
#' access to its public and private attributes and methods.
#'
#' @export
SampData_QA <- R6::R6Class(
  classname = "SampData_QA",
  lock_class = TRUE,
  inherit = SampData,
  public = list(
    #' @description
    #' Create a New `SampData_QA` Object
    #'
    #' @param params (`list`)
    #'   Parameters to initialize a [COGASData] object. Must also include a
    #'   `data.frame` or [tibble::tibble] named `cert_data` containing the
    #'   following columns:
    #'
    #'   * sample_id (`character`): sample ids
    #'   * c_cert (`numeric`): certified amount fractions
    #'   * uc_cert (`numeric`): uncertainties associated with the amount
    #'   fractions
    #'
    #' @return A new `SampData_QA` object.
    initialize = function(params) {
      # Expected columns to check `cert_data` against.
      expected_cols <- c(
        sample_id = "character",
        c_cert = "numeric",
        uc_cert = "numeric"
      )

      # Validating inputs.
      stopifnot(
        length(params) > 0,
        "cert_data" %in% names(params),
        is.data.frame(params$cert_data),
        purrr::imap_lgl(params$cert_data, \(col, name) {
          name %in% names(expected_cols) && class(col) == expected_cols[[name]]
        })
      )

      private$cert_data <- params$cert_data %>%
        dplyr::relocate(tidyselect::all_of(names(expected_cols))) %>%
        dplyr::arrange(c_cert) %>%
        dplyr::mutate(
          sample_id = forcats::as_factor(stringr::str_to_upper(sample_id))
        )

      super$initialize(params[names(params) != "cert_data"])
    },

    #' @description
    #' Print Object
    #'
    #' @return The calling object invisibly.
    print = function() {
      super$print()

      cat(
        "<SampData_QA>\n",
        "  Private Attributes:\n",
        "    certification data: tbl_df [",
        paste0(dim(private$cert_data), collapse = " x "),
        "    qa results plot: ",
        ifelse(is.null(private$qa_results_plot), "NULL\n", "ggplot\n"),
        "]\n",
        sep = ""
      )

      invisible(self)
    },

    #' @description
    #' Calculate Amount Fractions
    #'
    #' @details
    #' The parameters needed in the call to
    #' [`SampData$calc_amounts()`](SampData#method-calc_amounts) must be
    #' passed through `...`. The function additionally calculates differences
    #' between the certified amount fractions and the current calculations.
    #'
    #' @param ... Additional named arguments.
    #'
    #' @return The calling object invisibly.
    calc_amounts = function(...) {
      super$calc_amounts(...)

      private$amounts <- private$amounts %>%
        dplyr::right_join(private$cert_data, by = "sample_id") %>%
        dplyr::mutate(
          diff = amount - c_cert,
          udiff = sqrt(uamount^2 + uc_cert^2),
          percent_diff = 100 * diff / c_cert,
          zscore = diff / udiff
        )

      invisible(self)
    },

    #' @description
    #' Plot Data
    #'
    #' @details
    #' This function generates useful diagnostic plots for the data from the
    #' current analysis.
    #'
    #' @param units (`character(1)`)
    #'   Units to use for the y-axis label. Defaults to "mu*mol ~~ mol^-1".
    #'
    #' @return The calling object invisibly.
    plot_data = function(units = defaults$UNITS) {
      if (!is.null(private$amounts)) {
        super$plot_data(units)

        # Plotting amount fractions and expanded uncertainties of both the
        # certified and current values.
        private$qa_results_plot <- private$amounts %>%
          dplyr::mutate(uamount = 2 * uamount, uc_cert = 2 * uc_cert) %>%
          ggplot2::ggplot(ggplot2::aes(x = sample_id)) +
          ggplot2::geom_pointrange(
            ggplot2::aes(
              y = c_cert,
              ymin = c_cert - uc_cert,
              ymax = c_cert + uc_cert,
              shape = "Certified"
            ),
            size = 0.5
          ) +
          ggplot2::geom_pointrange(
            ggplot2::aes(
              y = amount,
              ymin = amount - uamount,
              ymax = amount + uamount,
              shape = "Current"
            ),
            size = 0.5
          ) +
          ggplot2::scale_shape_manual(
            name = NULL,
            values = c(1, 2),
            breaks = c("Certified", "Current")
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
    #' Get Data Frames
    #'
    #' @details
    #' This function returns a list of all private attributes that are data
    #' frames. If all attributes are not needed at the same time, it can be
    #' convenient to select a specific one from the function call, i.e.
    #' `get_data()$cleaned_data`.
    #'
    #' @return A `list` of attributes.
    get_data = function() {
      rlang::list2(!!!super$get_data(), cert_data = private$cert_data)
    },

    #' @description
    #' Get Plots
    #'
    #' @details
    #' This function returns a list of all private attributes that are plots.
    #' If all attributes are not needed at the same time, it can be convenient
    #' to select a specific one from the function call, i.e.
    #' `get_plots()$break_set_plot`.
    #'
    #' @return A `list` of attributes.
    get_plots = function() {
      rlang::list2(
        !!!super$get_plots(),
        qa_results_plot = private$qa_results_plot
      )
    }
  ),
  private = list(
    # [tibble::tibble] Certification data associated with analyzed samples.
    cert_data = NULL,

    # [ggplot2::ggplot] Plot of the certified and current amount fractions
    # for each sample.
    qa_results_plot = NULL
  )
)
