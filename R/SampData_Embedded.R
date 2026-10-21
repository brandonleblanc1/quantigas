#' [R6::R6Class] for Sample Data
#'
#' Derived class for processing data from COGAS. Specifically deals with data
#' for unknowns when the sample(s) of interest is(are) analyzed concurrently
#' with the calibration data. It inherits from the [SampData] class and thus
#' has access to its public and private attributes and methods.
#'
#' @export
SampData_Embedded <- R6::R6Class(
  classname = "SampData_Embedded",
  lock_class = TRUE,
  inherit = SampData,
  public = list(
    #' @description
    #' Clean Data
    #'
    #' @details
    #' This function reads in and cleans the data from the current analysis. It
    #' is explicitly reimplemented in this class to override the default
    #' [COGASData] method and specifically pull the data for `filt_samples` out
    #' of the calibration data.
    #'
    #' @param cal_obj ([CalData])
    #'   Object used for handling the calibration data.
    #' @param filt_samples (`character`)
    #'   Sample(s) of interest to filter from the calibration data.
    #' @param ... Unused arguments.
    #'
    #' @return The calling object invisibly.
    clean_data = function(cal_obj, filt_samples, ...) {
      cal_data <- cal_obj$get_data()

      private$raw_data <- cal_data$raw_data
      private$cleaned_data <- cal_data$cleaned_data %>%
        dplyr::filter(sample_id %in% filt_samples) %>%
        dplyr::mutate(dplyr::across(
          tidyselect::where(is.factor),
          forcats::fct_drop
        ))
      private$ratio_mc <- cal_obj$get_mc()$ratio_mc

      invisible(self)
    }
  )
)
