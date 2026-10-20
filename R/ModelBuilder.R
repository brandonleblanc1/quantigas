#' @title
#' ModelBuilder Class
#'
#' @description
#' Abstract base [R6::R6Class] for building different models for GSMG
#' measurement data. This class acts as a container for comparing multiple
#' candidate models, assessing the chosen model's performance, propagating
#' measurement uncertainties, and preparing visualizations of different types
#' of models.
#'
#' @examples
ModelBuilder <- R6::R6Class(
  classname = "ModelBuilder",
  lock_class = TRUE,
  private = list(
    #' @field cal_data (`data.frame()`)
    #' Calibration data being used to create the model.
    cal_data = NULL,

    #' @field models (`data.frame()`)
    #' Optimized parameters and statistical criteria for each candidate model.
    models = NULL,

    #' @description
    #' NLopt Wrapper
    #'
    #' @details
    #' This function is a wrapper for `nloptr()` to provide a common
    #' optimization interface for derived classes. Each derived class will
    #' implement its own `private$objective()` function corresponding to the
    #' specific model assumptions.
    #'
    #' @param theta0 (`numeric()`)
    #'   Initial values for the adjustable parameters of the model.
    #' @param lb (`numeric()`)
    #'   Lower bounds for the parameters.
    #' @param ub (`numeric()`)
    #'   Upper bounds for the parameters.
    #' @param func (`character(1)`)
    #'   Analysis function.
    #' @param x (`numeric()`)
    #'   Sampled x-values.
    #' @param y (`numeric()`)
    #'   Sampled y-values.
    #' @param ... Additional named arguments to pass to the `opts` parameter.
    #'
    #' @return `list()`
    optimizer = function(theta0, lb, ub, func, x, y, ...) {
      # Capturing additional options and giving default values for certain ones
      # if left unspecified.
      defaults <- list(
        algorithm = "NLOPT_LN_SBPLX",
        xtol_rel = 1e-7,
        maxtime = 60,
        maxeval = -1
      )

      nloptr::nloptr(
        x0 = theta0,
        eval_f = private$objective,
        lb = lb,
        ub = ub,
        opts = modifyList(defaults, list(...)),
        func = func,
        x = x,
        y = y
      )
    }
  ),
  public = list(
    #' @description
    #' Print [ModelBuilder] Object
    #'
    #' @return [ModelBuilder]
    print = function() {
      cat(
        "<ModelBuilder>\n",
        "  Private Attributes:\n",
        "    calibration data: tibble [",
        paste0(dim(private$cal_data), collapse = " x "),
        "]\n",
        "    candidate models data: tibble [",
        paste0(dim(private$models), collapse = " x "),
        "]\n",
        sep = ""
      )

      invisible(self)
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
    #' @return `list()`
    get_data = function() {
      list(cal_data = private$cal_data, models = private$models)
    },

    #' @description
    #' Save Object
    #'
    #' @details
    #' This function saves the object as an RDS file.
    #'
    #' @param filepath (`character(1)`)
    #'   Path to save the object to.
    #'
    #' @return [ModelBuilder]
    save_object = function(filepath) {
      saveRDS(self, file = filepath)

      invisible(self)
    }
  )
)

#' @title
#' Two-Standard Method (Linear Interpolation)
#'
#' @description
#' This function implements a linear interpolation method for predicting a
#' control sample (ratio of 1) based on the certified values of two standards.
#'
#' @param cal_obj ([CalData])
#'   Calibration data object.
#'
#' @return `data.frame()`
twostd_method <- function(cal_obj) {
  stopifnot(
    inherits(cal_obj, "CalData"),
    nrow(cal_obj$get_data()$cert_data) == 2
  )

  data <- cal_obj$get_data()

  cal_data <- data$cleaned_data %>%
    dplyr::right_join(data$cert_data, by = "sample_id") %>%
    dplyr::select(
      sample_id,
      sample_set,
      r = ratio,
      ur = uratio,
      c = c_cert,
      uc = uc_cert
    ) %>%
    dplyr::arrange(c, sample_set) %>%
    dplyr::mutate(
      dplyr::across(tidyselect::where(is.factor), forcats::fct_drop),
      sample_id = forcats::fct_reorder(sample_id, c),
      sample_set = forcats::fct_inorder(sample_set)
    )

  ratio_mc <- cal_obj$get_mc()$ratio_mc[
    as.numeric(levels(cal_data$sample_set)),
  ]
  amount_mc <- stats::rnorm(
    2 * ncol(ratio_mc),
    mean = cal_data$c[!duplicated(cal_data$sample_id)],
    sd = cal_data$uc[!duplicated(cal_data$sample_id)]
  ) %>%
    matrix(nrow = 2)

  ratio_splits <- split(seq_len(nrow(cal_data)), cal_data$sample_id)
  r1 <- t(ratio_mc[ratio_splits[[1]], ])
  r2 <- t(ratio_mc[ratio_splits[[2]], ])
  c1 <- amount_mc[1, ]
  c2 <- amount_mc[2, ]

  predicted_mc <- c1 + (1 - r1) * (c2 - c1) / (r2 - r1)

  tibble::tibble_row(
    sample_id = cal_obj$get_info()$control,
    amount = robustbase::huberM(predicted_mc)$mu,
    uamount = robustbase::Qn(predicted_mc),
    uamount_rel = 100 * uamount / amount
  )
}
