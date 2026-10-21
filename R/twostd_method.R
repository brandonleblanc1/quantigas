#' Two-Standard Method (Linear Interpolation)
#'
#' This function implements a linear interpolation method for predicting a
#' control sample (ratio of 1) based on the certified values of two standards.
#'
#' @param cal_obj ([CalData])
#'   Calibration data object.
#'
#' @return A [tibble::tibble] of results for the control sample.
#'
#' @export
twostd_method <- function(cal_obj) {
  stopifnot(
    inherits(cal_obj, "CalData"),
    nrow(cal_obj$get_data()$cert_data) == 2
  )

  data <- cal_obj$get_data()

  # Combining certification data with cleaned data.
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

  # Extracting portions of MC samples for each standard.
  ratio_splits <- split(seq_len(nrow(cal_data)), cal_data$sample_id)
  r1 <- t(ratio_mc[ratio_splits[[1]], ])
  r2 <- t(ratio_mc[ratio_splits[[2]], ])
  c1 <- amount_mc[1, ]
  c2 <- amount_mc[2, ]

  # Linear interpolation formula for predicting a value of 1. Note that R
  # performs vector operations on matrices in column-major order, which is why
  # the transpose of ratios is taken above. The vector then gets recycled across
  # the columns of the matrix.
  predicted_mc <- c1 + (1 - r1) * (c2 - c1) / (r2 - r1)

  # Summarizing results.
  tibble::tibble_row(
    sample_id = cal_obj$get_info()$control,
    amount = robustbase::huberM(predicted_mc)$mu,
    uamount = robustbase::Qn(predicted_mc),
    uamount_rel = 100 * uamount / amount
  )
}
