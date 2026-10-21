#' Calculate Gravimetric Results
#'
#' This function calculates differences between cylinders and a control sample
#' for gravimetric weighings when preparing PSMs. This may be made into a
#' separate [R6::R6Class] or [ModelBuilder] subclass in the future.
#'
#' @param filepath (`character(1)`)
#'   Path to the data file.
#' @param control (`character(1)`)
#'   Name of the control sample throughout the data file.
#'
#' @returns A [tibble::tibble] of results for the sample cylinders.
#'
#' @export
gravcalc <- function(filepath, control) {
  mass_obj <- COGASData$new(list(
    filepath = filepath,
    instrument = "BALANCE",
    compound = NA,
    control = control
  ))

  # Computing differences instead of ratios.
  diff_mc <- mass_obj$get_mc()$response_mc - mass_obj$get_mc()$control_mc

  cleaned_data <- mass_obj$get_data()$raw_data %>%
    dplyr::mutate(
      diff = response - control,
      udiff = apply(diff_mc, 1, robustbase::Qn)
    ) %>%
    dplyr::select(
      !c(tidyselect::ends_with("control"), tidyselect::ends_with("ratio"))
    ) %>%
    dplyr::filter_out(dplyr::when_any(
      sample_id == control,
      control_set == 0
    )) %>%
    dplyr::mutate(dplyr::across(
      tidyselect::where(is.factor),
      forcats::fct_drop
    )) %>%
    dplyr::relocate(tidyselect::ends_with("diff"), .after = day_set) %>%
    dplyr::arrange(sample_id, sample_set)

  # Summarizing results.
  diff_mc %>%
    tibble::as_tibble(
      .name_repair = "unique_quiet",
      rownames = "sample_set"
    ) %>%
    dplyr::left_join(cleaned_data, ., by = "sample_set") %>%
    tidyr::pivot_longer(tidyselect::starts_with("...")) %>%
    dplyr::summarise(
      diff = robustbase::huberM(value)$mu,
      udiff = robustbase::Qn(value),
      u_rel = 100 * udiff / diff,
      .by = sample_id
    )
}
