#' [R6::R6Class] for Model Building
#'
#' Base class for modelling GSMG measurement data. Essentially just a container
#' for comparing multiple candidate models, assessing model performance,
#' predicting unknown values, and preparing useful visualizations.
#'
#' @export
ModelBuilder <- R6::R6Class(
  classname = "ModelBuilder",
  lock_class = TRUE,
  public = list(
    #' @description
    #' Print useful information about the state of the object.
    #'
    #' @return The calling object invisibly.
    print = function() {
      cat(
        "<ModelBuilder>\n",
        "  Private Attributes:\n",
        "    calibration data: tbl_df [",
        paste0(dim(private$cal_data), collapse = " x "),
        "]\n",
        "    candidate models data: tbl_df [",
        paste0(dim(private$models), collapse = " x "),
        "]\n",
        sep = ""
      )

      invisible(self)
    },

    #' @description
    #' Get all private attributes that are data frames. If all attributes are
    #' not needed at the same time, it can be convenient to select a specific
    #' one from the function call, i.e. `get_data()$cal_data`.
    #'
    #' @return A `list` of attributes.
    get_data = function() {
      list(cal_data = private$cal_data, models = private$models)
    },

    #' @description
    #' Save the object to an RDS file.
    #'
    #' @param filepath (`character(1)`)
    #'   Path to save the object to.
    #'
    #' @return The calling object invisibly.
    save_object = function(filepath) {
      saveRDS(self, file = filepath)

      invisible(self)
    }
  ),
  private = list(
    # [tibble::tibble]: Calibration data being used to create the model.
    cal_data = NULL,

    # [tibble::tibble]: Optimized parameters and statistical criteria for
    # each candidate model.
    models = NULL,

    # This function is a wrapper for [nloptr::nloptr()] to provide a common
    # optimization interface for derived classes. Each derived class will
    # implement its own `private$objective()` function corresponding to the
    # specific model assumptions.
    optimizer = function(theta0, lb, ub, func, x, y, ...) {
      # Capturing additional options and giving default values for certain ones
      # if left unspecified.
      nlopts <- list(
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
        opts = modifyList(nlopts, list(...)),
        func = func,
        x = x,
        y = y
      )
    }
  )
)
