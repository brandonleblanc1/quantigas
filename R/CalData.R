#' [R6::R6Class] for Calibration Data
#'
#' Derived class for processing data from COGAS. Specifically deals with data
#' for instrument calibration and includes an attribute to store the
#' certification data associated with the calibration standards. Inherits from
#' the [COGASData] class and thus has access to its public and private
#' attributes and methods.
#'
#' @export
CalData <- R6::R6Class(
  classname = "CalData",
  lock_class = TRUE,
  inherit = COGASData,
  public = list(
    #' @description
    #' Create a new `CalData` object.
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
    #' @return A new `CalData` object.
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
    #' Print useful information about the state of the object.
    #'
    #' @return The calling object invisibly.
    print = function() {
      super$print()

      cat(
        "<CalData>\n",
        "  Private Attributes:\n",
        "    certification data: tbl_df [",
        paste0(dim(private$cert_data), collapse = " x "),
        "]\n",
        sep = ""
      )

      invisible(self)
    },

    #' @description
    #' Get all private attributes that are data frames. If all attributes are
    #' not needed at the same time, it can be convenient to select a specific
    #' one from the function call, i.e. `get_data()$cleaned_data`.
    #'
    #' @return A `list` of attributes.
    get_data = function() {
      rlang::list2(!!!super$get_data(), cert_data = private$cert_data)
    }
  ),
  private = list(
    # [tibble::tibble]: Certification data associated with calibration
    # standards.
    cert_data = NULL
  )
)
