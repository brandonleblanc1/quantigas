#' Two-Standard Method (Linear Interpolation)
#'
#' This function implements a linear interpolation method for predicting a
#' control sample (ratio of 1) based on the certified values of two standards.
#'
#' @param cal_obj ([CalData])
#'   Calibration data object.
#' @param filt_samples (`character`)
#'   Sample(s) of interest to exclude from the model. Defaults to an empty
#'   string.
#'
#' @return A [tibble::tibble] of results for the control sample.
#'
#' @export
twostd_method <- function(cal_obj, filt_samples = "") {
  stopifnot(
    inherits(cal_obj, "CalData"),
    length(setdiff(cal_obj$get_data()$cert_data$sample_id, filt_samples)) == 2
  )

  data <- cal_obj$get_data()

  # Combining certification data with cleaned data.
  cal_data <- data$cleaned_data %>%
    dplyr::right_join(data$cert_data, by = "sample_id") %>%
    dplyr::filter_out(sample_id %in% filt_samples) %>%
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
  ratio_splits <- split.data.frame(ratio_mc, rownames(ratio_mc))
  r1 <- t(ratio_splits[[1]])
  r2 <- t(ratio_splits[[2]])
  c1 <- amount_mc[1, ]
  c2 <- amount_mc[2, ]

  # Linear interpolation formula for predicting a value of 1. Note that R
  # performs vector operations on matrices in column-major order, which is why
  # the transpose of ratios is taken above. The vector then gets recycled across
  # the columns of the matrix.
  pred_mc <- c1 + (1 - r1) * (c2 - c1) / (r2 - r1)
  pred_sigma <- robustbase::covMcd(pred_mc)$cov

  # Calculating weighted average and corresponding variance.
  ones <- rep(1, nrow(pred_sigma))
  sigma_inv <- MASS::ginv(pred_sigma)
  denom <- as.numeric(t(ones) %*% sigma_inv %*% ones)
  w <- sigma_inv %*% ones / denom

  # Summarizing results.
  tibble::tibble_row(
    sample_id = cal_obj$get_info()$control,
    amount = as.numeric(t(w) %*% colMeans(pred_mc)),
    uamount = sqrt(1 / denom),
    uamount_rel = 100 * uamount / amount
  )
}
