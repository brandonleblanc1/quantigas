#' [R6::R6Class] for Model Building
#'
#' Derived class for modelling GSMG measurement data. Specifically deals with
#' building an error-in-variables model with all replicates data that includes
#' correlation within the ratios. Parameter uncertainties are obtained through
#' numerical methods. It inherits from the [PolyModel] class and thus has
#' access to its public and private attributes and methods.
#'
#' @export
PolyModel_Correlated <- R6::R6Class(
  classname = "PolyModel_Correlated",
  lock_class = TRUE,
  inherit = PolyModel,
  public = list(
    #' @description
    #' Create a new `PolyModel_Correlated` object.
    #'
    #' @param cal_obj ([CalData])
    #'   Calibration data object.
    #' @param filt_samples (`character`)
    #'   Sample(s) of interest to exclude from the model. Defaults to an empty
    #'   string.
    #'
    #' @return A new `PolyModel_Correlated` object.
    initialize = function(cal_obj, filt_samples = "") {
      stopifnot(inherits(cal_obj, "CalData"))

      cal_obj$summarize_data()
      cal_data <- cal_obj$get_data()
      cal_mc <- cal_obj$get_mc()

      # Combining certification data with cleaned data.
      private$cal_data <- cal_data$cleaned_data %>%
        dplyr::right_join(cal_data$cert_data, by = "sample_id") %>%
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

      # Getting initial values for rho and sigma, as well as lower and upper
      # bounds for rho.
      private$init_vals <- private$cal_data %>%
        dplyr::right_join(cal_data$summ_stats, by = "sample_id") %>%
        dplyr::select(
          sample_id,
          rho0 = avg_ratio,
          sigma0 = u_total,
          c,
          uc,
          n = n_repl
        ) %>%
        dplyr::summarise(
          dplyr::across(tidyselect::everything(), dplyr::last),
          rho_lb = rho0 - 6 * sigma0,
          rho_ub = rho0 + 6 * sigma0,
          .by = sample_id
        )

      private$n_std <- nrow(private$init_vals)
      private$n_repl <- sum(private$init_vals$n)
      private$nu <- private$n_repl - private$n_std

      # Subsetting ratio MC samples for relevant data.
      private$r_mc <- cal_mc$ratio_mc[
        as.numeric(levels(private$cal_data$sample_set)),
      ]

      private$r_sigma <- robustbase::covMcd(t(private$r_mc))$cov
    },

    #' @description
    #' Print useful information about the state of the object.
    #'
    #' @return The calling object invisibly.
    print = function() {
      super$print()

      cat(
        "<PolyModel_Replicates>\n",
        "  Private Attributes:\n",
        "    number standards: ",
        private$n_std,
        "\n",
        "    number replicates: ",
        private$n_repl,
        "\n",
        "    degrees of freedom: ",
        private$nu,
        "\n",
        "    selected model: ",
        private$selected_model,
        "\n",
        "    ratio covariances: matrix [",
        paste0(dim(private$r_sigma), collapse = " x "),
        "]\n",
        "    param covariances: matrix [",
        paste0(dim(private$p_sigma), collapse = " x "),
        "]\n",
        "    dark uncertainty mc samples: numeric [",
        length(private$lambda_mc),
        "]\n",
        "    initial values: tbl_df [",
        paste0(dim(private$init_vals), collapse = " x "),
        "]\n",
        sep = ""
      )

      invisible(self)
    },

    #' @description
    #' Build candidate polynomial models and run the optimization algorithm for
    #' each to assess their fit of the data. It also generates residual plots
    #' for each candidate model. The best model has the smallest signed value
    #' of AIC, BIC, or both; however, the model selected for use is usually the
    #' one in which the values of AIC and BIC show their first substantial drop
    #' as model complexity increases. These statistical criteria should be
    #' taken into account in conjunction with the residual plots, which should
    #' show no discernable patterns.
    #'
    #' @return The calling object invisibly.
    model_selection = function() {
      # Setting up for parallel processing.
      future::plan(future::multisession, workers = future::availableCores() - 2)

      private$models <- private$build_models() %>%
        dplyr::mutate(
          # Running the optimization algorithm for each candidate model.
          opt = furrr::future_map2(forms, funcs, \(form, func) {
            private$run_opt(form, func)
          })
        ) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          # Extracting parameters and calculating performance metrics.
          rho = list(opt$solution[1:private$n_std]),
          tau = exp(opt$solution[(private$n_std + 1)]),
          beta = list(opt$solution[-(1:(private$n_std + 1))]),
          xi = list(vapply(
            rho,
            \(r) eval(str2expression(funcs))(r, beta),
            numeric(1)
          )),
          r_resid = list(private$cal_data$r - rep(rho, private$init_vals$n)),
          r_zscore = list(r_resid / private$cal_data$ur),
          c_resid = list(private$init_vals$c - xi),
          c_zscore = list(c_resid / sqrt(private$init_vals$uc^2 + tau^2)),
          n_params = length(opt$solution),
          objective = opt$objective,
          aic = 2 * (n_params + opt$objective),
          bic = log(private$n_repl) * n_params + 2 * opt$objective,
          rmse = sqrt(mean(c_resid^2)),
          rse = sqrt(sum(c_resid^2) / (private$n_repl - n_params - 1)),
          r_sq = 1 -
            sum(c_resid^2) /
              sum((private$init_vals$c - mean(private$init_vals$c))^2),
          adj_r_sq = 1 -
            (1 - r_sq) * (private$n_repl - 1) / (private$n_repl - n_params - 1)
        ) %>%
        dplyr::ungroup()

      # Done with parallel processing.
      future::plan(future::sequential)

      # Generating plots of standardized residuals for each candidate model.
      private$r_resid_plot <- private$plot_models_resid("r")
      private$c_resid_plot <- private$plot_models_resid("c")

      invisible(self)
    },

    #' @description
    #' Perform a Monte Carlo simulation on the selected candidate model to
    #' propagate measurement uncertainty through the analysis function using
    #' the replicate MC samples of the calibration object. Assumes the certified
    #' amount fractions are normally distributed.
    #'
    #' @param model_num (`numeric(1)`)
    #'   Number identifying the candidate polynomial model to use in the
    #'   simulation. Defaults to 1.
    #' @param ... Additional options for the optimization routine.
    #'
    #' @return The calling object invisibly.
    simulate = function(model_num = defaults$MODEL_NUM, ...) {
      private$selected_model <- model_num

      # Extracting necessary parameters from the selected model.
      model_info <- private$models %>%
        dplyr::filter(.data$model_num == private$selected_model)

      solution <- model_info$opt[[1]]$solution
      K <- ncol(private$r_mc)

      # Approximating the second-partial derivatives of the objective function.
      hess <- numDeriv::hessian(
        \(x) {
          private$objective(
            theta = x,
            func = eval(str2expression(model_info$funcs)),
            x = private$cal_data$r,
            y = private$init_vals$c
          )
        },
        solution,
        method.args = list(eps = 1e-5, d = 1e-4, zero.tol = 1e-6, r = 10)
      )

      # Calculating the pseudo-inverse to get the parameter covariance matrix.
      # The pseudo-inverse is used instead of the direct inverse to handle
      # cases where strong correlation between parameters results in a singular
      # Hessian.
      private$p_sigma <- MASS::ginv(hess)

      # Drawing samples from the multivariate distribution defined by the
      # parameter covariance matrix.
      param_mc <- mvtnorm::rmvnorm(
        K,
        mean = solution,
        sigma = private$p_sigma
      ) %>%
        t()

      private$beta_mc <- param_mc[-(1:(private$n_std + 1)), ]

      # Generating random amount fraction samples from a normal distribution.
      private$c_mc <- stats::rnorm(
        private$n_std * K,
        mean = private$init_vals$c,
        sd = private$init_vals$uc
      ) %>%
        matrix(ncol = K, dimnames = list(private$init_vals$sample_id))

      # Generating random dark uncertainty samples from a normal distribution.
      private$lambda_mc <- stats::rnorm(K, sd = model_info$tau)

      subsets <- split(seq_len(private$n_repl), private$cal_data$sample_id)

      # Estimating residuals and uncertainties.
      r_resid <- do.call(
        rbind,
        purrr::map(seq_along(subsets), \(i) {
          sweep(private$r_mc[subsets[[i]], ], 2, param_mc[i, ])
        })
      )

      r_summ <- tibble::tibble(
        resid = private$init_vals$rho0 - model_info$rho[[1]],
        uresid = sqrt(
          private$init_vals$sigma0^2 + diag(private$p_sigma)[1:(private$n_std)]
        )
      )

      c_resid <- private$c_mc - self$predict(param_mc[1:(private$n_std), ])

      c_summ <- tibble::tibble(
        resid = private$init_vals$c - model_info$xi[[1]],
        uresid = sqrt(diag(robustbase::covMcd(t(c_resid))$cov))
      )

      # Summarizing residuals for both ratios and amounts fractions.
      private$sim_summ <- purrr::map(list(r = r_summ, c = c_summ), \(df) {
        df %>%
          dplyr::mutate(
            sample_id = private$init_vals$sample_id,
            zscore = resid / uresid
          ) %>%
          dplyr::relocate(sample_id)
      }) %>%
        dplyr::bind_rows(.id = "source") %>%
        tidyr::pivot_wider(
          names_from = source,
          values_from = c(resid, uresid, zscore),
          names_glue = "{source}_{.value}"
        ) %>%
        dplyr::relocate(
          tidyselect::starts_with("c"),
          .after = tidyselect::last_col()
        )

      # # Generating summary plots of residuals.
      private$r_mc_plot <- private$plot_mc_resid("r", r_resid)
      private$c_mc_plot <- private$plot_mc_resid("c", c_resid)

      invisible(self)
    },

    #' @description
    #' Predict unknowns utilizing the results obtained from
    #' [`simulate()`](#method-simulate). MC samples are returned and must be
    #' summarized by the caller.
    #'
    #' @param data (NULL | `numeric` | `matrix`)
    #'   Behavior will change based on the value provided:
    #'
    #'   * NULL: will calculate the value for a control sample by taking the
    #'   sums of each set of coefficients. Returns a vector.
    #'   * vector: will calculate the outer product of the vector with the
    #'   coefficients for each power of r. Returns a matrix.
    #'   * matrix: will calculate the product of each element of the matrix
    #'   with the corresponding coefficients for each power of r. Returns a
    #'   matrix.
    #'
    #' @return A `vector` or `matrix` with the calculated MC samples.
    predict = function(data = NULL) {
      if (is.null(private$beta_mc)) {
        self$simulate()
      }

      if (is.null(data)) {
        # Control sample is a ratio of 1, so just taking the sums of each pair
        # of coefficients.
        colSums(private$beta_mc) + private$lambda_mc
      } else {
        # Extracting the powers of r from the analysis function.
        powers <- eval(str2expression(stringr::str_extract(
          private$models$funcs[
            private$models$model_num == private$selected_model
          ],
          "c\\(.+\\)(?=\\))"
        )))

        seq_along(powers) %>%
          purrr::map(\(i) {
            if (is.matrix(data)) {
              sweep(data^powers[i], 2, private$beta_mc[i, ], FUN = "*")
            } else {
              outer(data^powers[i], private$beta_mc[i, ])
            }
          }) %>%
          purrr::reduce(`+`) %>%
          sweep(2, private$lambda_mc, FUN = "+")
      }
    },

    #' @description
    #' Get all private attributes that contain general information on the
    #' analysis. If all attributes are not needed at the same time, it can be
    #' convenient to select a specific one from the function call, i.e.
    #' `get_info()$n_std`.
    #'
    #' @return A `list` of attributes.
    get_info = function() {
      list(
        n_std = private$n_std,
        n_repl = private$n_repl,
        nu = private$nu,
        selected_model = private$selected_model
      )
    },

    #' @description
    #' Get all private attributes containing MC samples. If all attributes are
    #' not needed at the same time, it can be convenient to select a specific
    #' one from the function call, i.e. `get_mc()$r_mc`.
    #'
    #' @return A `list` of attributes.
    get_mc = function() {
      rlang::list2(
        !!!super$get_mc(),
        r_sigma = private$r_sigma,
        p_sigma = private$p_sigma,
        lambda_mc = private$lambda_mc
      )
    },

    #' @description
    #' Get all private attributes that are data frames. If all attributes are
    #' not needed at the same time, it can be convenient to select a specific
    #' one from the function call, i.e. `get_data()$cal_data`.
    #'
    #' @return A `list` of attributes.
    get_data = function() {
      rlang::list2(!!!super$get_data(), init_vals = private$init_vals)
    }
  ),
  private = list(
    # `numeric(1)`: Number of calibration standards.
    n_std = NULL,

    # `numeric(1)`: Total number of replicate response ratios.
    n_repl = NULL,

    # `numeric(1)`: Degrees of freedom.
    nu = NULL,

    # `numeric(1)`: Number identifying the model selected in the uncertainty
    # simulation.
    selected_model = NULL,

    # `matrix`: Covariance matrix of ratios.
    r_sigma = NULL,

    # `matrix`: Covariance matrix of parameters
    p_sigma = NULL,

    # `matrix`: Monte Carlo samples for the estimated dark uncertainty.
    lambda_mc = NULL,

    # [tibble::tibble]: Initial values to use for model setup.
    init_vals = NULL,

    # This function serves as the objective during optimization and represents
    # the likelihood of observing the data given the measurement model, the
    # calibration relationship, and the values of the model's adjustable
    # parameters.
    objective = function(theta, func, x, y) {
      # Extracting parameters from `theta`. `rho` is the true ratio means,
      # `tau` is the dark uncertainty lurking in the standards, and `beta` is
      # the polynomial coefficients.
      rho <- theta[1:private$n_std]
      tau <- exp(theta[(private$n_std + 1)])
      beta <- theta[-(1:(private$n_std + 1))]

      # Likelihood of observing the replicate ratios.
      loglik_r <- mvtnorm::dmvt(
        x - rep(rho, private$init_vals$n),
        sigma = private$r_sigma * (private$nu - 2) / private$nu,
        df = private$nu
      )

      xi <- vapply(rho, \(r) func(r, beta), numeric(1))
      uc_total <- sqrt(private$init_vals$uc^2 + tau^2)

      # Likelihood of observing the certified amount fractions.
      loglik_c <- stats::dnorm(y, mean = xi, sd = uc_total, log = TRUE)

      -sum(loglik_r, loglik_c)
    },

    # This function runs the optimization algorithm for the specified model.
    run_opt = function(form, func) {
      lm0 <- stats::lm(
        stats::formula(form),
        tibble::tibble(r = private$init_vals$rho0, c = private$init_vals$c)
      )
      beta0 <- unname(lm0$coefficients)
      tau0 <- stats::sigma(lm0)
      theta0 <- c(private$init_vals$rho0, log(tau0), beta0)

      lb <- c(private$init_vals$rho_lb, rep(-Inf, 1 + length(beta0)))
      ub <- c(private$init_vals$rho_ub, rep(Inf, 1 + length(beta0)))

      private$optimizer(
        theta0,
        lb,
        ub,
        eval(str2expression(func)),
        private$cal_data$r,
        private$init_vals$c,
        maxtime = 300
      )
    }
  )
)
