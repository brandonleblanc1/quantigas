#' [R6::R6Class] for Sample Data
#'
#' Derived class for processing data from COGAS. Specifically deals with PGVP
#' data for unknown samples and determining how accurate the vendor's value is.
#' It inherits from the [SampData] class and thus has access to its public and
#' private attributes and methods.
#'
#' @export
SampData_PGVP <- R6::R6Class(
  classname = "SampData_PGVP",
  lock_class = TRUE,
  inherit = SampData,
  public = list(
    #' @description
    #' Create a new `SampData_PGVP` object.
    #'
    #' @param params (`list`)
    #'   Parameters to initialize a [SampData] object. Must also include a
    #'   `data.frame` or [tibble::tibble] named `vendor_data` containing the
    #'   following columns:
    #'
    #'   * producer (`character`): producer names,
    #'   * pgvp_id (`character`): PGVP id for producer,
    #'   * sample_id (`character`): sample ids,
    #'   * vendor_amount (`numeric`): vendor amount fractions,
    #'   * vendor_u (`numeric`): expanded uncertainties associated with the
    #'   amount fractions.
    #'
    #' @return A new `SampData_PGVP` object.
    initialize = function(params) {
      # Expected columns to check `vendor_data` against.
      expected_cols <- c(
        producer = "character",
        pgvp_id = "character",
        sample_id = "character",
        vendor_amount = "numeric",
        vendor_u = "numeric"
      )

      # Validating inputs.
      stopifnot(
        length(params) > 0,
        "vendor_data" %in% names(params),
        is.data.frame(params$vendor_data),
        purrr::imap_lgl(params$vendor_data, \(col, name) {
          name %in% names(expected_cols) && class(col) == expected_cols[[name]]
        })
      )

      # Renaming the columns and setting private attributes.
      private$vendor_data <- params$vendor_data %>%
        dplyr::relocate(tidyselect::all_of(names(expected_cols))) %>%
        dplyr::arrange(pgvp_id) %>%
        dplyr::mutate(
          sample_id = stringr::str_to_upper(sample_id),
          dplyr::across(producer:sample_id, forcats::as_factor)
        )

      super$initialize(params[names(params) != "vendor_data"])
    },

    #' @description
    #' Print useful information about the state of the object.
    #'
    #' @return The calling object invisibly.
    print = function() {
      super$print()

      cat(
        "<SampData_PGVP>\n",
        "  Private Attributes:\n",
        "    vendor data: tbl_df [",
        paste0(dim(private$vendor_data), collapse = " x "),
        "]\n",
        "    vendor results: tbl_df [",
        paste0(dim(private$vendor_results), collapse = " x "),
        "]\n",
        "    percent diff plot: ",
        ifelse(is.null(private$percent_diff_plot), "NULL\n", "ggplot\n"),
        "    vendor results plot: ",
        ifelse(is.null(private$vendor_results_plot), "NULL\n", "ggplot\n"),
        sep = ""
      )

      invisible(self)
    },

    #' @description
    #' Calculate amount fractions and uncertainties for unknown samples, and
    #' additionally calculate percent differences of the vendor results to the
    #' calculated NIST results. The parameters needed in the call to
    #' [`SampData$calc_amounts()`](SampData#method-SampData-calc_amounts)
    #' must be passed through `...`.
    #'
    #' @param ... Additional named arguments.
    #'
    #' @return The calling object invisibly.
    calc_amounts = function(...) {
      super$calc_amounts(...)

      private$vendor_results <- private$amounts %>%
        dplyr::select(!uamount_rel) %>%
        dplyr::right_join(private$vendor_data, ., by = "sample_id") %>%
        dplyr::relocate(day_set, .after = sample_id) %>%
        dplyr::mutate(
          uamount = 2 * uamount,
          percent_diff = 100 *
            dplyr::if_else(
              is.na(producer),
              (amount - vendor_amount) / vendor_amount,
              (vendor_amount - amount) / amount
            )
        )

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
      if (!is.null(private$vendor_results)) {
        super$plot_data(units)

        # Plotting percent differences between the vendor assigned value and the
        # NIST calculated value by sample id.
        private$percent_diff_plot <- private$vendor_results %>%
          ggplot2::ggplot(ggplot2::aes(x = sample_id, y = percent_diff)) +
          ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
          ggplot2::geom_hline(
            ggplot2::aes(yintercept = 2.2),
            color = "red",
            linetype = 2
          ) +
          ggplot2::geom_hline(
            ggplot2::aes(yintercept = -2.2),
            color = "red",
            linetype = 2
          ) +
          ggplot2::geom_point(size = 2) +
          ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2)) +
          ggplot2::labs(x = NULL, y = "Percent Difference") +
          ggplot2::theme_bw()

        # Plotting amount fractions and expanded uncertainties of both the
        # vendor results and NIST results.
        private$vendor_results_plot <- private$vendor_results %>%
          ggplot2::ggplot(ggplot2::aes(x = sample_id)) +
          ggplot2::geom_pointrange(
            ggplot2::aes(
              y = vendor_amount,
              ymin = vendor_amount - vendor_u,
              ymax = vendor_amount + vendor_u,
              shape = "Vendor"
            ),
            size = 0.5
          ) +
          ggplot2::geom_pointrange(
            ggplot2::aes(
              y = amount,
              ymin = amount - uamount,
              ymax = amount + uamount,
              shape = "NIST"
            ),
            size = 0.5
          ) +
          ggplot2::scale_shape_manual(
            name = NULL,
            values = c(1, 2),
            breaks = c("NIST", "Vendor")
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
    #' Get all private attributes that are data frames. If all attributes are
    #' not needed at the same time, it can be convenient to select a specific
    #' one from the function call, i.e. `get_data()$cleaned_data`.
    #'
    #' @return A `list` of attributes.
    get_data = function() {
      rlang::list2(
        !!!super$get_data(),
        vendor_data = private$vendor_data,
        vendor_results = private$vendor_results
      )
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
        percent_diff_plot = private$percent_diff_plot,
        vendor_results_plot = private$vendor_results_plot
      )
    }
  ),
  private = list(
    # [tibble::tibble]: Vendor data for cylinders being analyzed.
    vendor_data = NULL,

    # [tibble::tibble]: Results of vendor analyses compared to NIST analyses.
    vendor_results = NULL,

    # [ggplot2::ggplot]: Plot of the calculated percent differences for each
    # sample.
    percent_diff_plot = NULL,

    # [ggplot2::ggplot]: Plot of the vendor and NIST amount fractions for
    # each sample.
    vendor_results_plot = NULL
  )
)
